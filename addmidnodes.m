function [tri6,xnod,ynod]=addmidnodes(tri3,xv,yv)
% Add edge midpoint nodes to a triangular mesh, giving P2 connectivity.
% tri6 columns: 3 vertices, then midpoints of edges (1,2),(2,3),(3,1).
nv=length(xv); nele=size(tri3,1);
ed=[tri3(:,[1 2]);tri3(:,[2 3]);tri3(:,[3 1])];
[ue,~,ic]=unique(sort(ed,2),'rows');
xnod=[xv;(xv(ue(:,1))+xv(ue(:,2)))/2];
ynod=[yv;(yv(ue(:,1))+yv(ue(:,2)))/2];
tri6=[tri3,nv+reshape(ic,nele,3)];
end
