function [l2,h1s,bl2] = errnorms2(tri6,tri3,xnod,ynod,p,d)
% P2 version of errnorms: L2(Omega), H1-seminorm, L2(contact boundary)
% norms of the nodal field d.
enod = [1,2;2,3;3,1];
nele = size(tri6,1);
l2 = 0; h1s = 0; bl2 = 0;
for iel = 1:nele
    iv = tri6(iel,:)'; xc = xnod(tri3(iel,:)); yc = ynod(tri3(iel,:));
    [gcx,gcy,gv] = trigauc(xc,yc,5);
    for j = 1:length(gv)
        [fi,fix,fiy,area] = basescalar2(gcx(j),gcy(j),xc,yc);
        l2 = l2 + gv(j)*area*(fi'*d(iv))^2;
        h1s = h1s + gv(j)*area*((fix'*d(iv))^2 + (fiy'*d(iv))^2);
    end
end
for im = 1:size(p,2)
    iel = p(im).iel; inei = p(im).inei;
    iv = tri6(iel,:); xc = xnod(tri3(iel,:)); yc = ynod(tri3(iel,:));
    nod = tri3(iel,enod(inei,:));
    xcut = xnod(nod); ycut = ynod(nod);
    dl = sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
    ngau = 4; [gc,gw,~] = gauss(ngau,0,1);
    for j = 1:ngau
        xm = gc(j)*xcut(1)+(1-gc(j))*xcut(2);
        ym = gc(j)*ycut(1)+(1-gc(j))*ycut(2);
        [fi,~,~,~] = basescalar2(xm,ym,xc,yc);
        bl2 = bl2 + gw(j)*dl*(fi'*d(iv))^2;
    end
end
l2 = sqrt(l2); h1s = sqrt(h1s); bl2 = sqrt(bl2);
end
