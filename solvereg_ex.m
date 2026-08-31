function [u,iite,reshist] = solvereg_ex(tri,tri3,xnod,ynod,p,deg,pow,gamma0,exfun,ind)
% Newton solver for the regularized Nitsche contact method on a general
% mesh, with obstacle g = 0 on the contact boundary described by p,
% inhomogeneous Dirichlet data taken from the exact solution on the nodes
% ind, load f from the exact solution, and regularization epsilon = dl^pow.
% deg = 1 (P1) or 2 (P2); tri holds the element dofs, tri3 the vertices.
% NOTE: volume integrals use the element area (not detj = 2*area as in the
% older routines), so gamma0 here has its nominal meaning.
enod=[1,2;2,3;3,1];
nele=size(tri,1); nno=length(xnod); nmult=size(p,2); neq=nno;
ndof=size(tri,2); nsize=nele*ndof^2;
qvol=5; ngau=4;
u=zeros(nno,1);
[uD,~,~,~]=exfun(xnod(ind),ynod(ind));
u(ind)=uD;                       % Dirichlet data, kept fixed by the Newton loop
int=setdiff(1:neq,ind);
F=sparse(neq,1); F1=F; mite=200; reshist=zeros(mite,1);
nstep=0;

% --- a_h : volume stiffness and load ---
[row,col,val,up]=assemble(nsize);
for iel=1:nele
    iv=tri(iel,:)'; xc=xnod(tri3(iel,:)); yc=ynod(tri3(iel,:));
    [gcx,gcy,gv]=trigauc(xc,yc,qvol);
    sele=zeros(ndof); fele=zeros(ndof,1);
    for ng=1:length(gv)
        [fi,fix,fiy,area]=basisfun(deg,gcx(ng),gcy(ng),xc,yc);
        [~,~,~,fq]=exfun(gcx(ng),gcy(ng));
        bw=[fix';fiy'];
        fele=fele+gv(ng)*area*fi*fq;
        sele=sele+gv(ng)*area*(bw)'*bw;
    end
    [row,col,val,up]=assemble(sele,iv,row,col,val,up);
    F(iv)=F(iv)+fele;
end
% --- a_h : -1/gamma (sigma_n,sigma_n) on the contact boundary ---
for im=1:nmult
    iel=p(im).iel; inei=p(im).inei;
    iv=tri(iel,:); xc=xnod(tri3(iel,:)); yc=ynod(tri3(iel,:));
    nod=tri3(iel,enod(inei,:));
    xcut=xnod(nod); ycut=ynod(nod);
    dl=sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
    gamma=gamma0/dl; nvec=[ycut(2)-ycut(1),xcut(1)-xcut(2)]/dl;
    [gc,gw,~]=gauss(ngau,0,1);
    sele=zeros(ndof);
    for j=1:ngau
        xm=gc(j)*xcut(1)+(1-gc(j))*xcut(2);
        ym=gc(j)*ycut(1)+(1-gc(j))*ycut(2);
        [~,fix,fiy,~]=basisfun(deg,xm,ym,xc,yc);
        fin=nvec(1)*fix+nvec(2)*fiy;
        sele=sele-(1/gamma)*gw(j)*dl*(fin)*fin';
    end
    [row,col,val,up]=assemble(sele,iv(:),row,col,val,up);
end
ahmat=sparse(row(1:up),col(1:up),val(1:up),neq,neq);

% --- Newton loop ---
for iite=1:mite
    [row,col,val,up]=assemble(nsize);
    for im=1:nmult
        iel=p(im).iel; inei=p(im).inei;
        iv=tri(iel,:); xc=xnod(tri3(iel,:)); yc=ynod(tri3(iel,:));
        nod=tri3(iel,enod(inei,:));
        xcut=xnod(nod); ycut=ynod(nod);
        dl=sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
        gamma=gamma0/dl; ginv=1/gamma;
        nvec=[ycut(2)-ycut(1),xcut(1)-xcut(2)]/dl;
        [gc,gw,~]=gauss(ngau,0,1);
        sele=zeros(ndof); fele=zeros(ndof,1);
        epsilon=dl^pow;                    % epsilon = 4s
        sqrt_epsilon=sqrt(epsilon);
        for j=1:ngau
            xm=gc(j)*xcut(1)+(1-gc(j))*xcut(2);
            ym=gc(j)*ycut(1)+(1-gc(j))*ycut(2);
            [fi,fix,fiy,~]=basisfun(deg,xm,ym,xc,yc);
            fin=nvec(1)*fix+nvec(2)*fiy;
            ui=dot(fi,u(iv)); uin=dot(fin,u(iv));
            test=fi-ginv*fin;
            w=ui-ginv*uin;
            pp=hypot(w,sqrt_epsilon);
            if w>=0
                phi=(w+pp)/2;
            else
                phi=epsilon/(2*(pp-w));
            end
            Dphi=phi/pp;
            fele=fele-gamma*gw(j)*dl*phi*test;
            sele=sele+gamma*gw(j)*dl*Dphi*test*(test');
        end
        F1(iv)=F1(iv)+fele;
        [row,col,val,up]=assemble(sele,iv(:),row,col,val,up);
    end
    Smat=sparse(row(1:up),col(1:up),val(1:up),neq,neq);
    Smat=Smat+ahmat; f=F+F1; F1=0*F1;
    res=f-ahmat*u;
    resn=norm(res(int));
    reshist(iite)=resn;                          % residual before the update
    if iite==1, res0=resn; end
    if resn<max(1e-11,1e-12*res0), break; end   % relative safeguard on fine meshes
    v=full(Smat(int,int)\res(int));
    du=zeros(neq,1); du(int)=v;
    u=u+du;
    nstep=nstep+1;
end
neval=iite;
reshist=reshist(1:neval);
iite=nstep;
end
