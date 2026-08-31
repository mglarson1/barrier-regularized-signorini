function [u,iite,resn] = solvereg_eps(nodes,xnod,ynod,p,pow)
% Newton solver for the regularized Nitsche contact method (as solvereg.m),
% with regularization epsilon = dl^pow (dl = local edge length).
ind = find(xnod == min(xnod));
gamma0 = 10;
enod = [1,2;2,3;3,1];
nele = size(nodes,1); nno = length(xnod); nmult = size(p,2); neq = nno;
u = zeros(nno,1);
nsize = nele*3*3;
F = sparse(neq,1); F1 = F;
mite = 200;
nstep = 0;

% --- a_h : volume stiffness + load ---
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
% --- a_h : -1/gamma (sigma_n, sigma_n) on contact edges ---
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
        [~,fix,fiy,~] = basescalar(xm,ym,xc,yc);
        fin = nvec(1)*fix+nvec(2)*fiy;
        sele = sele - (1/gamma)*gw(j)*dl*(fin)*fin';
    end
    [row,col,val,up] = assemble(sele,iv(:),row,col,val,up);
end
ahmat = sparse(row(1:up),col(1:up),val(1:up),neq,neq);

int = setdiff(1:neq,ind);
% --- Newton loop ---
for iite = 1:mite
    [row,col,val,up] = assemble(nsize);
    for im = 1:nmult
        iel = p(im).iel; inei = p(im).inei;
        iv = nodes(iel,:); xc = xnod(iv); yc = ynod(iv);
        nod = iv(enod(inei,:));
        xcut = xnod(nod); ycut = ynod(nod);
        dl = sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
        gamma = gamma0/dl; ginv = 1/gamma;
        nvec = [ycut(2)-ycut(1),xcut(1)-xcut(2)]/dl;
        ngau = 2; [gc,gw,~] = gauss(ngau,0,1);
        sele = zeros(3); fele = zeros(3,1);
        epsilon = dl^pow;                    % epsilon = 4s
        sqrt_epsilon = sqrt(epsilon);
        for j = 1:ngau
            xm = gc(j)*xcut(1)+(1-gc(j))*xcut(2);
            ym = gc(j)*ycut(1)+(1-gc(j))*ycut(2);
            wei = gw(j);
            [fi,fix,fiy,~] = basescalar(xm,ym,xc,yc);
            fin = nvec(1)*fix+nvec(2)*fiy;
            ui = dot(fi,u(iv)); uin = dot(fin,u(iv));
            test = fi - ginv*fin;
            w = ui - ginv*uin;
            pp = hypot(w,sqrt_epsilon);
            if w >= 0
                phi = (w + pp)/2;
            else
                phi = epsilon/(2*(pp - w));
            end
            Dphi = phi/pp;                    % varphi_{+,s}' = varphi_{+,s}/pp
            fele = fele - gamma*wei*dl*phi*test;
            sele = sele + gamma*wei*dl*Dphi*test*(test');
        end
        F1(iv) = F1(iv) + fele;
        [row,col,val,up] = assemble(sele,iv(:),row,col,val,up);
    end
    Smat = sparse(row(1:up),col(1:up),val(1:up),neq,neq);
    Smat = Smat + ahmat; f = F + F1; F1 = 0*F1;
    res = f - ahmat*u;
    resn = norm(res(int));
    if iite == 1, res0 = resn; end
    if resn < max(1e-10,1e-12*res0), break; end   % relative safeguard on fine meshes
    v = full(Smat(int,int)\res(int));
    du = zeros(neq,1); du(int) = v;
    u = u + du;
    nstep = nstep + 1;
end
iite = nstep;
end
