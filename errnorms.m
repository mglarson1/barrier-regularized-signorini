function [l2,h1s,bl2] = errnorms(nodes,xnod,ynod,p,d)
% L2(Omega), H1-seminorm, and L2(contact boundary) norms of nodal field d
enod = [1,2;2,3;3,1];
nele = size(nodes,1);
l2 = 0; h1s = 0; bl2 = 0;
for iel = 1:nele
    iv = nodes(iel,:)'; xc = xnod(iv); yc = ynod(iv);
    [gcx,gcy,gv] = trigauc(xc,yc,2);
    for j = 1:length(gv)
        [fi,fix,fiy,area] = basescalar(gcx(j),gcy(j),xc,yc);
        l2 = l2 + gv(j)*area*(fi'*d(iv))^2;
    end
    h1s = h1s + area*((fix'*d(iv))^2 + (fiy'*d(iv))^2);
end
for im = 1:size(p,2)
    iel = p(im).iel; inei = p(im).inei;
    iv = nodes(iel,:); xc = xnod(iv); yc = ynod(iv);
    nod = iv(enod(inei,:));
    xcut = xnod(nod); ycut = ynod(nod);
    dl = sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
    ngau = 2; [gc,gw,~] = gauss(ngau,0,1);
    for j = 1:ngau
        xm = gc(j)*xcut(1)+(1-gc(j))*xcut(2);
        ym = gc(j)*ycut(1)+(1-gc(j))*ycut(2);
        [fi,~,~,~] = basescalar(xm,ym,xc,yc);
        bl2 = bl2 + gw(j)*dl*(fi'*d(iv))^2;
    end
end
l2 = sqrt(l2); h1s = sqrt(h1s); bl2 = sqrt(bl2);
end
