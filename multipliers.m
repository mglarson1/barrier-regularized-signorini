function p=multipliers(trisq,neighbors,xnodsq,ynodsq)
nele=size(trisq,1);
nod=[1,2;2,3;3,1];
i1=1;
for iel=1:nele
    iv=trisq(iel,:);
    for in=1:3
        xn=xnodsq(iv(nod(in,:)));
        yn=ynodsq(iv(nod(in,:)));
        if(neighbors(iel,in)==0)
            if(max(xn)>min(xnodsq))
            p(i1).iel=iel;p(i1).inei=in;i1=i1+1;
            end
        end
    end
end

            
            
