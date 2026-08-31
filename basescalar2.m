function [fi,fix,fiy,area]=basescalar2(x,y,xc,yc)
% P2 (quadratic) basis functions on a triangle.
% xc,yc = coordinates of the 3 VERTICES. Local node order:
% 1,2,3 = vertices; 4 = mid(1,2), 5 = mid(2,3), 6 = mid(3,1).
[L,Lx,Ly,area]=basescalar(x,y,xc,yc);
fi =[L(1)*(2*L(1)-1); L(2)*(2*L(2)-1); L(3)*(2*L(3)-1); ...
     4*L(1)*L(2); 4*L(2)*L(3); 4*L(3)*L(1)];
fix=[(4*L(1)-1)*Lx(1); (4*L(2)-1)*Lx(2); (4*L(3)-1)*Lx(3); ...
     4*(Lx(1)*L(2)+L(1)*Lx(2)); 4*(Lx(2)*L(3)+L(2)*Lx(3)); 4*(Lx(3)*L(1)+L(3)*Lx(1))];
fiy=[(4*L(1)-1)*Ly(1); (4*L(2)-1)*Ly(2); (4*L(3)-1)*Ly(3); ...
     4*(Ly(1)*L(2)+L(1)*Ly(2)); 4*(Ly(2)*L(3)+L(2)*Ly(3)); 4*(Ly(3)*L(1)+L(3)*Ly(1))];
end
