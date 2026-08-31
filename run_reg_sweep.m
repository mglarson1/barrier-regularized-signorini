% Convergence of the regularization error as a function of the exponent in
% s ~ h^pow. Theory (Section 5, central path argument):
%   ||u_s - u_0||_H1 <~ (gamma*s)^(1/2) = O(h^((pow-1)/2)),
% and optimal order requires rate >= k = 1, i.e. pow >= 2k+1 = 3 for P1.
% (The older direct estimate predicted rate pow/2 - 1 and required pow >= 4.)
levels = [8 16 32 64 128 256];
pows   = [2 3 4 5];
H1 = zeros(numel(levels),numel(pows));
IT = zeros(numel(levels),numel(pows));
IT0 = zeros(numel(levels),1);
for li = 1:numel(levels)
    n = levels(li);
    [tri,xnod,ynod] = unitsquaremesh(n);
    TR = triangulation(tri,xnod,ynod); nb = TR.neighbors; nb(isnan(nb)) = 0;
    nb = [nb(:,3),nb(:,1),nb(:,2)];
    p = multipliers(tri,nb,xnod,ynod);
    [u0,it0] = solvefixed_q(tri,xnod,ynod,p);    % unregularized, s = 0
    IT0(li) = it0;
    for pi = 1:numel(pows)
        [us,it] = solvereg_eps(tri,xnod,ynod,p,pows(pi));
        [~,h1,~] = errnorms(tri,xnod,ynod,p,us-u0);
        H1(li,pi) = h1; IT(li,pi) = it;
    end
end
fprintf('\nRegularization error ||u_s - u_0||_H1,  s ~ h^pow/4\n');
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
fprintf('\nNewton iterations\n%6s %4s','n','it0'); fprintf('  s~h^%-1d',pows); fprintf('\n');
for li = 1:numel(levels)
    fprintf('%6d %4d',levels(li),IT0(li)); fprintf('  %5d',IT(li,:)); fprintf('\n');
end
save('results_p1.mat','H1','IT','IT0','levels','pows');
