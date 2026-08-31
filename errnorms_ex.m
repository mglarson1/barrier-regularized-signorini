function [l2,h1] = errnorms_ex(tri,tri3,xnod,ynod,deg,uh,exfun)
% L2 and full H1 norm of the error u - uh against the exact solution.
qvol=5; l2=0; h1s=0;
for iel=1:size(tri,1)
    iv=tri(iel,:)'; xc=xnod(tri3(iel,:)); yc=ynod(tri3(iel,:));
    [gcx,gcy,gv]=trigauc(xc,yc,qvol);
    for j=1:length(gv)
        [fi,fix,fiy,area]=basisfun(deg,gcx(j),gcy(j),xc,yc);
        [ue,uex,uey,~]=exfun(gcx(j),gcy(j));
        l2 =l2 +gv(j)*area*(ue-fi'*uh(iv))^2;
        h1s=h1s+gv(j)*area*((uex-fix'*uh(iv))^2+(uey-fiy'*uh(iv))^2);
    end
end
l2=sqrt(l2); h1=sqrt(l2^2+h1s);
end
