function [u,iite,resn] = solvefixed_q(nodes,xnod,ynod,p)
% Quiet version of solvefixed.m: unregularized (s=0) semismooth/active-set
% reference solution.
ind = find(xnod == min(xnod));
BC = zeros(size(ind));
gamma0 = 10;
enod = [1,2;2,3;3,1];
nele = size(nodes,1); nno = length(xnod); nmult = size(p,2); neq = nno;
u = zeros(nno,1); u(ind) = BC;
nsize = nele*3*3;
F = sparse(neq,1); F1 = F;
mite = 200;
nstep = 0;

[row,col,val,up] = assemble(nsize);
for iel = 1:nele
    iv = nodes(iel,:)'; xc = xnod(iv); yc = ynod(iv);
    [gcx,gcy,gv] = trigauc(xc,yc,1);
    sele = zeros(3); fele = zeros(3,1);
    for ng = 1:length(gv)
        [fi,fix,fiy,detj] = basescalar(gcx(ng),gcy(ng),xc,yc);
        bweps = [fix';fiy'];
        fele = fele + gv(ng)*detj*fi*floa(gcx(ng),gcy(ng));
        sele = sele + gv(ng)*detj*(bweps)'*bweps;
    end
    [row,col,val,up] = assemble(sele,iv,row,col,val,up);
    F(iv) = F(iv) + fele;
end
Gmat = sparse(row(1:up),col(1:up),val(1:up),neq,neq);

int = setdiff(1:neq,ind);
for iite = 1:mite
    [row,col,val,up] = assemble(nsize);
    for im = 1:nmult
        iel = p(im).iel; inei = p(im).inei;
        iv = nodes(iel,:); xc = xnod(iv); yc = ynod(iv);
        nod = iv(enod(inei,:));
        xcut = xnod(nod); ycut = ynod(nod);
        dl = sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
        gamma = gamma0/dl;
        nvec = [ycut(2)-ycut(1),xcut(1)-xcut(2)]/dl;
        ngau = 2; [gc,gw,~] = gauss(ngau,0,1);
        sele = zeros(3);
        for j = 1:ngau
            xm = gc(j)*xcut(1)+(1-gc(j))*xcut(2);
            ym = gc(j)*ycut(1)+(1-gc(j))*ycut(2);
            wei = gw(j);
            [fi,fix,fiy,~] = basescalar(xm,ym,xc,yc);
            fin = nvec(1)*fix+nvec(2)*fiy;
            ui = dot(fi,u(iv)); uin = dot(fin,u(iv));
            sigaug = double(ui - (1/gamma)*uin > 0);   % Heaviside(P(u))
            sele = sele - (1/gamma)*wei*dl*(fin)*fin';
            sele = sele + sigaug*wei*dl*gamma*(fi-(1/gamma)*fin)*(fi-(1/gamma)*fin)';
        end
        [row,col,val,up] = assemble(sele,iv(:),row,col,val,up);
    end
    Smat = sparse(row(1:up),col(1:up),val(1:up),neq,neq);
    Smat = Smat + Gmat; f = F + F1;
    res = f - Smat*u;
    resn = norm(res(int));
    if iite == 1, res0 = resn; end
    if resn < max(1e-12,1e-12*res0), break; end
    v = full(Smat(int,int)\f(int));
    du = zeros(neq,1); du(int) = v;
    uold = u; u = du;
    nstep = nstep + 1;
    if norm(u-uold) <= 1e-14*max(1,norm(u)), break; end   % active set unchanged
end
iite = nstep;
end
