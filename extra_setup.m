function [tri,tri3,xn,yn,p,ind] = extra_setup(n,deg)
% Mesh, contact edges and Dirichlet nodes for the exact-solution examples on
% (-1,1)x(0,1), with the contact boundary the bottom edge y = 0 and Dirichlet
% data on the remaining three sides.  Same construction as run_exact_all.m.
[tri3,xv,yv]=rectmesh(2*n,n,[-1 1],[0 1]);
if deg==1
    tri=tri3; xn=xv; yn=yv;
else
    [tri,xn,yn]=addmidnodes(tri3,xv,yv);
end
TR=triangulation(tri3,xv,yv); nb=TR.neighbors; nb(isnan(nb))=0;
nb=[nb(:,3),nb(:,1),nb(:,2)];
p=multipliers_bottom(tri3,nb,xn,yn,0);
ind=find(abs(xn+1)<1e-10|abs(xn-1)<1e-10|abs(yn-1)<1e-10);
end
