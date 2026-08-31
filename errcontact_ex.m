function [pl2,viol,aset,gapc,penl2,boundl2] = errcontact_ex(tri,tri3,xnod,ynod,p,deg,pow,gamma0,uh,exfun)
% Quantities on the contact boundary for the regularized Nitsche method,
% needed for the remarks on the contact pressure and on strict feasibility.
%
%   pl2   L2(Gamma_C) norm of sigma_n(u) - lambda_h, where
%         lambda_h = -gamma*phi_{+,s}(P(u_h)) is the discrete contact
%         pressure.  Remark 5.5 bounds ||sigma_n(u_mu) - lambda_h|| by
%         gamma0^{1/2} h^{k-1/2}; here the reference is the *exact* pressure
%         sigma_n(u), so the measured quantity also carries the distance
%         from the central path pressure to sigma_n(u), which the remark
%         does not estimate.
%   viol  max over the contact quadrature points of the constraint
%         violation max(u_h - g, 0), with g = 0.  Remark 5.3 predicts this
%         to be of the order of the finite element error.
%   aset  measure of the computed active set {P(u_h) > 0}.  The exact
%         contact set is [-1,0] in both examples, of measure 1.
%   gapc  largest artificial separation g - u_h over the quadrature points
%         that lie in the exact contact set, identified by sigma_n(u) < 0.
%         Remark 5.3 predicts a separation mu/|sigma_n(u_mu)| when s is
%         chosen too large.
%   penl2 L2(Gamma_C) norm of the penetration (u_h-g)_+.
%   boundl2 L2(Gamma_C) norm of gamma^{-1}|sigma_n(u_h)-lambda_h|.
%         The pointwise algebraic inequality implies penl2 <= boundl2 and
%         therefore provides a like-norm verification of the feasibility
%         estimate without comparing a maximum norm with an L2 norm.
%
% The quadrature, the normal and the regularization follow solvereg_ex
% exactly, so that lambda_h here is the same quantity the solver assembles.
% Note s = epsilon/4 with epsilon = dl^pow, since the solver evaluates
% phi = (w + sqrt(w^2 + epsilon))/2 = w/2 + sqrt(w^2/4 + epsilon/4).
enod=[1,2;2,3;3,1]; ngau=4;
pl2=0; viol=0; aset=0; gapc=0; penl2=0; boundl2=0;
[gc,gw,~]=gauss(ngau,0,1);
for im=1:size(p,2)
    iel=p(im).iel; inei=p(im).inei;
    iv=tri(iel,:); xc=xnod(tri3(iel,:)); yc=ynod(tri3(iel,:));
    nod=tri3(iel,enod(inei,:));
    xcut=xnod(nod); ycut=ynod(nod);
    dl=sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
    gamma=gamma0/dl; ginv=1/gamma;
    nvec=[ycut(2)-ycut(1),xcut(1)-xcut(2)]/dl;
    epsilon=dl^pow;                         % epsilon = 4s
    sqrt_epsilon=sqrt(epsilon);
    for j=1:ngau
        xm=gc(j)*xcut(1)+(1-gc(j))*xcut(2);
        ym=gc(j)*ycut(1)+(1-gc(j))*ycut(2);
        [fi,fix,fiy,~]=basisfun(deg,xm,ym,xc,yc);
        fin=nvec(1)*fix+nvec(2)*fiy;
        ui=dot(fi,uh(iv)); uin=dot(fin,uh(iv));
        w=ui-ginv*uin;                       % P(u_h)
        q=hypot(w,sqrt_epsilon);
        if w>=0
            phi=(w+q)/2;
        else
            phi=epsilon/(2*(q-w));
        end
        lam=-gamma*phi;                      % discrete contact pressure
        [~,uex,uey,~]=exfun(xm,ym);
        sn=nvec(1)*uex+nvec(2)*uey;          % sigma_n(u), same normal
        pl2=pl2+gw(j)*dl*(sn-lam)^2;
        pen=max(ui,0);                       % g = 0
        dbound=abs(uin-lam)/gamma;
        penl2=penl2+gw(j)*dl*pen^2;
        boundl2=boundl2+gw(j)*dl*dbound^2;
        viol=max(viol,ui);
        if w>0, aset=aset+gw(j)*dl; end
        if sn<-1e-12, gapc=max(gapc,-ui); end
    end
end
pl2=sqrt(pl2); penl2=sqrt(penl2); boundl2=sqrt(boundl2);
viol=max(viol,0);
end
