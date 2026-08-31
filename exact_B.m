function [u,ux,uy,f]=exact_B(x,y)
% Classical Signorini solution with generic free boundary behaviour, g = 0,
% contact boundary y = 0:
%   u = -r^{3/2} cos(3*theta/2),  f = 0  (harmonic)
% On y = 0: u = 0 and d_n u = -(3/2)sqrt(r) < 0 for x < 0 (contact),
%           u = -|x|^{3/2} < 0 and d_n u = 0 for x > 0 (free).
% Free boundary at the origin; u is in H^{5/2-eps} only.
r=sqrt(x.^2+y.^2); th=atan2(y,x);
u  = -r.^1.5.*cos(1.5*th);
ux = -1.5*sqrt(r).*cos(0.5*th);
uy =  1.5*sqrt(r).*sin(0.5*th);
f  = zeros(size(x));
end
