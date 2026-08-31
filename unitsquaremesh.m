function [tri,xnod,ynod] = unitsquaremesh(n)
% n x n uniform quad grid split into triangles as in maintri.m
[X,Y] = meshgrid(linspace(0,1,n+1));
xnod = X(:); ynod = Y(:);
id = @(i,j) (j-1)*(n+1) + i;   % column-major: i = y-index, j = x-index
tri = zeros(2*n*n,3); k = 0;
for j = 1:n
    for i = 1:n
        n1 = id(i,j); n2 = id(i,j+1); n3 = id(i+1,j+1); n4 = id(i+1,j);
        tri(k+1,:) = [n1,n2,n3];
        tri(k+2,:) = [n1,n3,n4];
        k = k + 2;
    end
end
end
