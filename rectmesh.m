function [tri,xnod,ynod]=rectmesh(nx,ny,xlim,ylim)
% Uniform triangulation of the rectangle xlim x ylim with nx by ny cells,
% each split into two triangles. Elements are counterclockwise (detj > 0).
[X,Y]=meshgrid(linspace(xlim(1),xlim(2),nx+1),linspace(ylim(1),ylim(2),ny+1));
xnod=X(:); ynod=Y(:);
id=@(i,j) (j-1)*(ny+1)+i;      % i = y-index, j = x-index (column major)
tri=zeros(2*nx*ny,3); k=0;
for j=1:nx
    for i=1:ny
        n1=id(i,j); n2=id(i,j+1); n3=id(i+1,j+1); n4=id(i+1,j);
        tri(k+1,:)=[n1,n2,n3];
        tri(k+2,:)=[n1,n3,n4];
        k=k+2;
    end
end
end
