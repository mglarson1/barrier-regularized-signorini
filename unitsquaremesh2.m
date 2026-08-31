function [tri6,tri3,xnod,ynod]=unitsquaremesh2(n)
% P2 mesh on the unit square: vertex mesh from unitsquaremesh plus edge
% midpoint nodes. tri6 columns: 3 vertices, then midpoints of edges
% (1,2), (2,3), (3,1).
[tri3,xv,yv]=unitsquaremesh(n);
nv=length(xv);
nele=size(tri3,1);
ed=[tri3(:,[1 2]);tri3(:,[2 3]);tri3(:,[3 1])];
[ue,~,ic]=unique(sort(ed,2),'rows');
xm=(xv(ue(:,1))+xv(ue(:,2)))/2;
ym=(yv(ue(:,1))+yv(ue(:,2)))/2;
xnod=[xv;xm]; ynod=[yv;ym];
mid=nv+reshape(ic,nele,3);
tri6=[tri3,mid];
end
