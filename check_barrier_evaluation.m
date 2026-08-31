% Verify the stable evaluation of the smoothed positive part over a range
% that includes the smallest regularization parameters used in the paper.
% The defining identities are
%   phi_s(w) = (w + sqrt(w^2+4s))/2,
%   phi_s(w)*(phi_s(w)-w) = s,
%   phi_s'(w) = phi_s(w)/sqrt(w^2+4s).

wvals=[-1,-1e-3,-1e-12,0,1e-12,1e-3,1];
svals=10.^(-(1:4:29));
max_comp_error=0;
max_derivative_error=0;

for is=1:numel(svals)
    s=svals(is);
    epsilon=4*s;
    sqrt_epsilon=sqrt(epsilon);
    for iw=1:numel(wvals)
        w=wvals(iw);
        q=hypot(w,sqrt_epsilon);
        if w>=0
            phi=(w+q)/2;
            dual_slack=epsilon/(2*(q+w));
        else
            phi=epsilon/(2*(q-w));
            dual_slack=(q-w)/2;
        end
        dphi=phi/q;

        comp_error=abs(phi*dual_slack-s)/s;
        derivative_error=abs(dphi+dual_slack/q-1);
        max_comp_error=max(max_comp_error,comp_error);
        max_derivative_error=max(max_derivative_error,derivative_error);
    end
end

assert(max_comp_error<5e-14,'Complementarity identity lost accuracy.');
assert(max_derivative_error<5e-14,'Derivative identity lost accuracy.');
fprintf('Stable barrier evaluation verified: complementarity %.3e, derivative %.3e.\n', ...
    max_comp_error,max_derivative_error);
