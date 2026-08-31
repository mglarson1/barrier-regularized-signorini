% Additional numerical evidence for Section 6 of the scalar paper:
%   (1) Newton residual histories and the observed order of convergence,
%       plus a robustness sweep in s with full steps (Remark 4.1);
%   (2) convergence of the discrete contact pressure (Remark 5.5);
%   (4) the violation of the constraint (Remark 5.3);
%   (5) the artificial gap and the shrinking active set when s is too
%       large (Remark 5.3).
% Output: numerics_extra.mat, read by make_extra_figs.m
%
% Geometry, elements and parameters follow run_exact_all.m: the domain is
% (-1,1)x(0,1), the contact boundary the bottom edge, gamma0 = 10*deg, and
% the regularization is s = h^pow/4.

EX={'A',@exact_A;'B',@exact_B};

% ---------- (1) Newton residual histories, fixed mesh ----------
% At the threshold s = h^{2k+1}/4 and well below it, to expose the
% shrinking radius of fast local convergence predicted by Remark 4.1.
RES=struct();
for deg=[1 2]
    n=32; gamma0=10*deg; k=deg;
    [tri,tri3,xn,yn,p,ind]=extra_setup(n,deg);
    pows=[2*k+1, 2*k+3, 2*k+5, 2*k+7];
    H=cell(1,numel(pows));
    for ip=1:numel(pows)
        [~,it,rh]=solvereg_ex(tri,tri3,xn,yn,p,deg,pows(ip),gamma0,@exact_A,ind);
        H{ip}=rh;
        fprintf('residuals P%d  s = h^%2d/4 : %2d iterations\n',deg,pows(ip),it);
    end
    RES.(sprintf('p%d',deg))=struct('n',n,'pows',pows,'hist',{H});
end

% ---------- (1b) robustness in s with full steps ----------
% The Newton system stays positive definite for every s > 0, but the radius
% of quadratic convergence shrinks.  Full steps are taken throughout, so a
% failure here is a failure of the undamped iteration, not of solvability.
ROB=struct();
for deg=[1 2]
    n=32; gamma0=10*deg; k=deg;
    [tri,tri3,xn,yn,p,ind]=extra_setup(n,deg);
    pows=(2*k+1):2:(2*k+15);
    IT=nan(size(pows)); OK=false(size(pows));
    for ip=1:numel(pows)
        [~,it,rh]=solvereg_ex(tri,tri3,xn,yn,p,deg,pows(ip),gamma0,@exact_A,ind);
        IT(ip)=it; OK(ip)= rh(end) < max(1e-11,1e-12*rh(1));
        fprintf('robust  P%d  s = h^%2d/4 : it %3d  converged %d\n',deg,pows(ip),it,OK(ip));
    end
    ROB.(sprintf('p%d',deg))=struct('n',n,'pows',pows,'IT',IT,'OK',OK);
end

% ---------- (2)+(4) contact pressure and constraint violation ----------
levels={[8 16 32 64 128],[8 16 32 64]};
CON=struct();
for ic=1:2
    for deg=[1 2]
        k=deg; pow=2*k+1; gamma0=10*deg; lv=levels{deg};
        PL=zeros(size(lv)); VI=PL; AS=PL; PN=PL; BD=PL;
        for li=1:numel(lv)
            n=lv(li);
            [tri,tri3,xn,yn,p,ind]=extra_setup(n,deg);
            uh=solvereg_ex(tri,tri3,xn,yn,p,deg,pow,gamma0,EX{ic,2},ind);
            [PL(li),VI(li),AS(li),~,PN(li),BD(li)]= ...
                errcontact_ex(tri,tri3,xn,yn,p,deg,pow,gamma0,uh,EX{ic,2});
        end
        CON.(sprintf('%s%d',EX{ic,1},deg))=struct( ...
            'h',1./lv,'pl2',PL,'viol',VI,'aset',AS,'penl2',PN,'boundl2',BD);
        fprintf('contact %s P%d done\n',EX{ic,1},deg);
    end
end

% ---------- (5) artificial gap and active set versus s ----------
% Fixed mesh, s swept from far above the threshold down to far below it.
% For large s the central path is strictly feasible by a wide margin and the
% body lifts off the obstacle: the active set shrinks and the gap grows.
GAP=struct();
for deg=[1 2]
    n=32; gamma0=10*deg; k=deg;
    [tri,tri3,xn,yn,p,ind]=extra_setup(n,deg);
    pows=1:2:(2*k+11);                  % increasing exponent = decreasing s
    AS=zeros(size(pows)); GC=AS; VI=AS; PN=AS; BD=AS; IT=AS;
    for ip=1:numel(pows)
        [uh,it]=solvereg_ex(tri,tri3,xn,yn,p,deg,pows(ip),gamma0,@exact_A,ind);
        [~,VI(ip),AS(ip),GC(ip),PN(ip),BD(ip)]= ...
            errcontact_ex(tri,tri3,xn,yn,p,deg,pows(ip),gamma0,uh,@exact_A);
        IT(ip)=it;
        fprintf('gap P%d  s = h^%2d/4 : active %.4f  gap %.3e  it %d\n', ...
            deg,pows(ip),AS(ip),GC(ip),it);
    end
    GAP.(sprintf('p%d',deg))=struct('n',n,'pows',pows,'s',(1/n).^pows/4, ...
        'aset',AS,'gapc',GC,'viol',VI,'penl2',PN,'boundl2',BD,'IT',IT);
end

save('numerics_extra.mat','RES','ROB','CON','GAP');
disp('numerics_extra.mat saved');
