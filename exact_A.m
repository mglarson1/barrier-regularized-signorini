function [u,ux,uy,f]=exact_A(x,y)
% Manufactured exact solution of the boundary obstacle problem, g = 0,
% contact boundary y = 0:
%   u = -x_+^3 + y*(-x)_+^3,   f = -Laplace(u)
% On y = 0: u = -x_+^3 <= 0, d_n u = -(-x)_+^3 <= 0, and u*d_n u = 0.
% Contact set x <= 0, free boundary at x = 0. u is C^{2,1}, hence in H^3.
xp=max(x,0); xm=max(-x,0);
u  = -xp.^3 + y.*xm.^3;
ux = -3*xp.^2 - 3*y.*xm.^2;
uy = xm.^3;
f  = 6*xp - 6*y.*xm;
end
