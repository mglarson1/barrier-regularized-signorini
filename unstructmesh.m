function [tri,xnod,ynod] = unstructmesh(n,seed)
% Unstructured triangulation of (-1,1)x(0,1) with nominal element size 1/n,
% for the exact-solution examples.  Boundary nodes are placed at spacing 1/n
% so that the contact edge y = 0 and the three Dirichlet sides are resolved
% exactly; the interior nodes are a randomly jittered grid.  The domain is
% convex, so the Delaunay triangulation fills it exactly.  A fixed seed makes
% the mesh reproducible.  Elements are returned counterclockwise, as in
% rectmesh.
if nargin<2, seed=0; end
rng(seed);
h=1/n;
xb=linspace(-1,1,2*n+1)'; yb=linspace(0,1,n+1)';
P=[ xb, zeros(size(xb));            % bottom, y = 0, the contact boundary
    xb, ones(size(xb));             % top
   -ones(size(yb)), yb;             % left
    ones(size(yb)), yb ];           % right
P=unique(round(P,12),'rows');
[gx,gy]=meshgrid(-1+h:h:1-h, h:h:1-h);
gx=gx(:); gy=gy(:); j=0.35*h;
gx=gx+j*(2*rand(size(gx))-1);
gy=gy+j*(2*rand(size(gy))-1);
keep = gx>-1+0.3*h & gx<1-0.3*h & gy>0.3*h & gy<1-0.3*h;
P=[P; gx(keep), gy(keep)];
DT=delaunayTriangulation(P);
tri=DT.ConnectivityList; xnod=DT.Points(:,1); ynod=DT.Points(:,2);
a=(xnod(tri(:,2))-xnod(tri(:,1))).*(ynod(tri(:,3))-ynod(tri(:,1))) ...
 -(xnod(tri(:,3))-xnod(tri(:,1))).*(ynod(tri(:,2))-ynod(tri(:,1)));
fl=a<0; tri(fl,[2 3])=tri(fl,[3 2]);
end
