function p=multipliers_bottom(tri3,neighbors,xnod,ynod,y0,tol)
% Contact edges = boundary edges lying on the line y = y0.
if nargin<6, tol=1e-10; end
nod=[1,2;2,3;3,1]; i1=0; p=struct('iel',{},'inei',{});
for iel=1:size(tri3,1)
    iv=tri3(iel,:);
    for in=1:3
        if neighbors(iel,in)==0
            yn=ynod(iv(nod(in,:)));
            if all(abs(yn-y0)<tol)
                i1=i1+1; p(i1).iel=iel; p(i1).inei=in;
            end
        end
    end
end
end
