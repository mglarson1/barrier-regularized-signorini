% P2 version of run_reg_sweep: convergence of the regularization error for
% s ~ h^pow/4. Theory (Section 5, central path argument):
%   ||u_s - u_0||_H1 <~ (gamma*s)^(1/2) = O(h^((pow-1)/2)),
% and optimal order requires rate >= k = 2, i.e. pow >= 2k+1 = 5.
% (The older direct estimate required pow >= 2(k+1) = 6.)
levels = [8 16 32 64 128];
pows   = [3 4 5 6 7];
H1 = zeros(numel(levels),numel(pows));
IT = zeros(numel(levels),numel(pows));
IT0 = zeros(numel(levels),1);
for li = 1:numel(levels)
    n = levels(li);
    [tri6,tri3,xnod,ynod] = unitsquaremesh2(n);
    TR = triangulation(tri3,xnod(1:max(tri3(:))),ynod(1:max(tri3(:))));
    nb = TR.neighbors; nb(isnan(nb)) = 0;
    nb = [nb(:,3),nb(:,1),nb(:,2)];
    p = multipliers(tri3,nb,xnod,ynod);
    [u0,it0,res0] = solvefixed_q2(tri6,tri3,xnod,ynod,p);
    IT0(li) = it0;
    for pi = 1:numel(pows)
        [us,it] = solvereg_eps2(tri6,tri3,xnod,ynod,p,pows(pi));
        [~,h1,~] = errnorms2(tri6,tri3,xnod,ynod,p,us-u0);
        H1(li,pi) = h1; IT(li,pi) = it;
    end
end
fprintf('\nP2: regularization error ||u_s - u_0||_H1,  s ~ h^pow/4\n');
fprintf('%6s','n'); fprintf('   s~h^%-1d  EOC ',pows); fprintf('\n');
for li = 1:numel(levels)
    fprintf('%6d',levels(li));
    for pi = 1:numel(pows)
        if li == 1, e = NaN; else, e = log2(H1(li-1,pi)/H1(li,pi)); end
        fprintf(' %8.2e %4.2f',H1(li,pi),e);
    end
    fprintf('\n');
end
fprintf('%6s','pred'); fprintf('          %4.2f',(pows-1)/2); fprintf('   <- predicted rate (pow-1)/2\n');
fprintf('\nNewton iterations (it0 = reference)\n%6s %4s','n','it0'); fprintf('  s~h^%-1d',pows); fprintf('\n');
for li = 1:numel(levels)
    fprintf('%6d %4d',levels(li),IT0(li)); fprintf('  %5d',IT(li,:)); fprintf('\n');
end
save('results_p2.mat','H1','IT','IT0','levels','pows');
