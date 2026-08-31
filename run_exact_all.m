% Regenerate and store all convergence data for the figures of Section 6:
% examples A and B, P1 and P2, plus the effect of the regularization parameter.
levels={[8 16 32 64 128 256],[8 16 32 64 128]};
cases={'A',@exact_A;'B',@exact_B};
DATA=struct();
for ic=1:2
  for deg=[1 2]
    k=deg; pow=2*k+1; gamma0=10*deg; lv=levels{deg};
    L2=zeros(size(lv)); H1=L2; IT=L2;
    for li=1:numel(lv)
        n=lv(li);
        [tri3,xv,yv]=rectmesh(2*n,n,[-1 1],[0 1]);
        if deg==1, tri=tri3; xn=xv; yn=yv; else, [tri,xn,yn]=addmidnodes(tri3,xv,yv); end
        TR=triangulation(tri3,xv,yv); nb=TR.neighbors; nb(isnan(nb))=0;
        nb=[nb(:,3),nb(:,1),nb(:,2)];
        p=multipliers_bottom(tri3,nb,xn,yn,0);
        ind=find(abs(xn+1)<1e-10|abs(xn-1)<1e-10|abs(yn-1)<1e-10);
        [uh,it]=solvereg_ex(tri,tri3,xn,yn,p,deg,pow,gamma0,cases{ic,2},ind);
        [L2(li),H1(li)]=errnorms_ex(tri,tri3,xn,yn,deg,uh,cases{ic,2});
        IT(li)=it;
    end
    DATA.(sprintf('%s%d',cases{ic,1},deg))=struct('h',1./lv,'L2',L2,'H1',H1,'IT',IT);
    fprintf('%s P%d done\n',cases{ic,1},deg);
  end
end
% effect of the regularization parameter, P1 and example A
lv=[8 16 32 64]; SV=struct();
for pow=[3 5 7]
    L2=zeros(size(lv)); H1=L2;
    for li=1:numel(lv)
        n=lv(li);
        [tri3,xv,yv]=rectmesh(2*n,n,[-1 1],[0 1]);
        TR=triangulation(tri3,xv,yv); nb=TR.neighbors; nb(isnan(nb))=0;
        nb=[nb(:,3),nb(:,1),nb(:,2)];
        p=multipliers_bottom(tri3,nb,xv,yv,0);
        ind=find(abs(xv+1)<1e-10|abs(xv-1)<1e-10|abs(yv-1)<1e-10);
        uh=solvereg_ex(tri3,tri3,xv,yv,p,1,pow,10,@exact_A,ind);
        [L2(li),H1(li)]=errnorms_ex(tri3,tri3,xv,yv,1,uh,@exact_A);
    end
    SV.(sprintf('p%d',pow))=struct('h',1./lv,'L2',L2,'H1',H1);
    fprintf('s ~ h^%d done\n',pow);
end
save('exact_conv_data.mat','DATA','SV');
disp('saved');
