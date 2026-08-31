% Two further studies for Section 6 of the scalar paper:
%   (6) the Nitsche parameter gamma0, whose admissible range the theorem
%       leaves as "sufficiently large", gamma0 > C_I;
%   (7) unstructured meshes, on which the constant C_I of the inverse
%       inequality is sensitive to element quality.
% Output: numerics_gamma_mesh.mat, read by make_extra_figs.m
%
% Both use Example A at the threshold s = h^{2k+1}/4.

% ---------- (6) sweep in gamma0 ----------
g0=[0.5 1 1.5 2 3 4 5 7 10 15 20 30 50 100 200 500];
GAM=struct();
for deg=[1 2]
    k=deg; pow=2*k+1; lv=[16 32 64];
    H1=nan(numel(lv),numel(g0)); IT=H1; OK=false(numel(lv),numel(g0));
    for li=1:numel(lv)
        n=lv(li);
        [tri,tri3,xn,yn,p,ind]=extra_setup(n,deg);
        for ig=1:numel(g0)
            try
                [uh,it,rh]=solvereg_ex(tri,tri3,xn,yn,p,deg,pow,g0(ig),@exact_A,ind);
                if all(isfinite(uh))
                    [~,H1(li,ig)]=errnorms_ex(tri,tri3,xn,yn,deg,uh,@exact_A);
                    IT(li,ig)=it;
                    OK(li,ig)= rh(end) < max(1e-11,1e-12*rh(1));
                end
            catch
                % singular or indefinite system: leave the entry as NaN
            end
        end
        fprintf('gamma0 sweep P%d  n=%3d done\n',deg,n);
    end
    GAM.(sprintf('p%d',deg))=struct('g0',g0,'lv',lv,'H1',H1,'IT',IT,'OK',OK);
end

% ---------- (7) unstructured meshes ----------
% Same quantities as run_exact_all.m, on Delaunay meshes of jittered points.
% hmax is the largest element diameter, reported alongside the nominal 1/n.
EX={'A',@exact_A;'B',@exact_B};
levels={[8 16 32 64],[8 16 32]};
UNS=struct();
for ic=1:2
    for deg=[1 2]
        k=deg; pow=2*k+1; gamma0=10*deg; lv=levels{deg};
        L2=zeros(size(lv)); H1=L2; IT=L2; HM=L2; PL=L2; VI=L2;
        for li=1:numel(lv)
            n=lv(li);
            [tri3,xv,yv]=unstructmesh(n);
            if deg==1, tri=tri3; xn=xv; yn=yv; else, [tri,xn,yn]=addmidnodes(tri3,xv,yv); end
            TR=triangulation(tri3,xv,yv); nb=TR.neighbors; nb(isnan(nb))=0;
            nb=[nb(:,3),nb(:,1),nb(:,2)];
            p=multipliers_bottom(tri3,nb,xn,yn,0);
            ind=find(abs(xn+1)<1e-10|abs(xn-1)<1e-10|abs(yn-1)<1e-10);
            e=[tri3(:,[1 2]);tri3(:,[2 3]);tri3(:,[3 1])];
            HM(li)=max(sqrt((xv(e(:,1))-xv(e(:,2))).^2+(yv(e(:,1))-yv(e(:,2))).^2));
            [uh,it]=solvereg_ex(tri,tri3,xn,yn,p,deg,pow,gamma0,EX{ic,2},ind);
            [L2(li),H1(li)]=errnorms_ex(tri,tri3,xn,yn,deg,uh,EX{ic,2});
            [PL(li),VI(li)]=errcontact_ex(tri,tri3,xn,yn,p,deg,pow,gamma0,uh,EX{ic,2});
            IT(li)=it;
            fprintf('unstructured %s P%d  n=%3d  hmax=%.4f  H1=%.3e  it=%d\n', ...
                EX{ic,1},deg,n,HM(li),H1(li),it);
        end
        UNS.(sprintf('%s%d',EX{ic,1},deg))= ...
            struct('h',1./lv,'hmax',HM,'L2',L2,'H1',H1,'IT',IT,'pl2',PL,'viol',VI);
    end
end

save('numerics_gamma_mesh.mat','GAM','UNS');
disp('numerics_gamma_mesh.mat saved');
