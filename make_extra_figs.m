% Figures for the additional numerical evidence of Section 6.
% Input: numerics_extra.mat (run_numerics_extra.m), results_p1/p2.mat (run_reg_sweep*.m)
% Output: figures/newton-convergence.pdf   (iteration counts + residual histories)
%         figures/pressure-convergence.pdf
%         figures/gap-sweep.pdf
set_publication_graphics;
S=load('numerics_extra.mat');
S1=load('results_p1.mat'); S2=load('results_p2.mat');
fnt='Times New Roman';
fig=[fullfile(fileparts(mfilename('fullpath')),'figures') filesep];
if ~isfolder(fig), mkdir(fig); end
c=lines(6); mk={'o','s','^','d','v'};

% ---------- Newton: iteration counts (top) and residual histories (bottom) ----------
f=figure('Units','centimeters','Position',[2 2 18 13],'Visible','off');
for d=1:2
    subplot(2,2,d)
    if d==1, Q=S1; else, Q=S2; end
    for k=1:numel(Q.pows)
        semilogx(Q.levels,Q.IT(:,k),'-o','Color',c(k,:),'LineWidth',1.2, ...
            'MarkerSize',5,'MarkerFaceColor',c(k,:)); hold on
    end
    if isfield(Q,'IT0')
        semilogx(Q.levels,Q.IT0,'k--s','LineWidth',1.2,'MarkerSize',5);
    end
    grid on, box on, set(gca,'FontSize',9,'FontName',fnt)
    xlabel('\itn\rm = 1/\ith','FontSize',10,'FontName',fnt)
    if d==1, ylabel('Newton iterations','FontSize',10,'FontName',fnt); end
    ylim([6 24]); xlim([7 300])       % headroom so the legend clears the curves
    lg=arrayfun(@(a) sprintf('\\its\\rm = \\ith\\rm^{%d}/4',a),Q.pows,'UniformOutput',false);
    if isfield(Q,'IT0'), lg=[lg,{'active set, \its\rm = 0'}]; end
    set(legend(lg,'Location','best','FontSize',7,'FontName',fnt),'ItemTokenSize',[12 6])
    title(sprintf('\\itP\\rm_{%d}',d),'FontSize',10,'FontName',fnt,'FontWeight','normal')

    subplot(2,2,2+d)
    R=S.RES.(sprintf('p%d',d));
    lg=cell(1,numel(R.pows));
    for ip=1:numel(R.pows)
        r=R.hist{ip};
        semilogy(0:numel(r)-1,r,['-' mk{ip}],'Color',c(ip,:),'LineWidth',1.2, ...
            'MarkerSize',4.5,'MarkerFaceColor',c(ip,:)); hold on
        lg{ip}=sprintf('\\its\\rm = \\ith\\rm^{%d}/4',R.pows(ip));
    end
    grid on, box on, set(gca,'FontSize',9,'FontName',fnt)
    xlabel('Newton step','FontSize',10,'FontName',fnt)
    if d==1, ylabel('residual','FontSize',10,'FontName',fnt); end
    ylim([1e-16 1e2]); xlim([0 20])
    set(legend(lg,'Location','best','FontSize',7,'FontName',fnt),'ItemTokenSize',[12 6])
    title(sprintf('\\itP\\rm_{%d},  \\itn\\rm = %d',d,R.n), ...
        'FontSize',10,'FontName',fnt,'FontWeight','normal')
end
exportgraphics(f,[fig 'newton-convergence.pdf'],'ContentType','vector','BackgroundColor','white'); close(f)

% ---------- contact pressure convergence ----------
f=figure('Units','centimeters','Position',[2 2 18 7.5],'Visible','off');
for panel=1:2
    subplot(1,2,panel); k=panel;
    A=S.CON.(sprintf('A%d',k)); B=S.CON.(sprintf('B%d',k));
    loglog(A.h,A.pl2,'-o','Color',c(1,:),'LineWidth',1.3,'MarkerSize',5, ...
        'MarkerFaceColor',c(1,:)); hold on
    loglog(B.h,B.pl2,'-s','Color',c(2,:),'LineWidth',1.3,'MarkerSize',5, ...
        'MarkerFaceColor',c(2,:));
    h=A.h;
    loglog(h,A.pl2(end)*(h/h(end)).^(k-0.5),'k--','LineWidth',1.0);
    loglog(h,A.pl2(end)*(h/h(end)).^k,'k-','LineWidth',1.0);
    grid on, box on, set(gca,'FontSize',9,'FontName',fnt)
    xlabel('\ith','FontSize',10,'FontName',fnt)
    if panel==1
        ylabel('\it||\sigma_n(u) - \lambda_h||\rm_{\Gamma_C}', ...
            'FontSize',10,'FontName',fnt);
    end
    yl=ylim; ylim([yl(1) yl(2)*6]);      % headroom for the legend
    set(legend({'A','B',sprintf('\\ith\\rm^{%g}  (central path)',k-0.5), ...
        sprintf('\\ith\\rm^{%d}  (observed)',k)}, ...
        'Location','northwest','FontSize',7.5,'FontName',fnt),'ItemTokenSize',[12 6])
    title(sprintf('\\itP\\rm_{%d},  \\its\\rm = \\ith\\rm^{%d}/4',k,2*k+1), ...
        'FontSize',10,'FontName',fnt,'FontWeight','normal')
end
exportgraphics(f,[fig 'pressure-convergence.pdf'],'ContentType','vector','BackgroundColor','white'); close(f)

% ---------- artificial gap: active set (left) and gap (right) ----------
f=figure('Units','centimeters','Position',[2 2 18 7.5],'Visible','off');
G1=S.GAP.p1; G2=S.GAP.p2;
sthr=[(1/G1.n)^3/4, (1/G2.n)^5/4];        % threshold s = h^{2k+1}/4
for panel=1:2
    subplot(1,2,panel)
    for d=1:2
        G=S.GAP.(sprintf('p%d',d));
        if panel==1, y=G.aset; else, y=G.gapc; end
        semilogx(G.s,y,['-' mk{d}],'Color',c(d,:),'LineWidth',1.3, ...
            'MarkerSize',5,'MarkerFaceColor',c(d,:)); hold on
    end
    for d=1:2
        xline(sthr(d),'--','Color',c(d,:),'LineWidth',1.0);
    end
    grid on, box on, set(gca,'FontSize',9,'FontName',fnt)
    xlabel('\its','FontSize',10,'FontName',fnt)
    if panel==1
        plot(xlim,[1 1],':','Color',[.4 .4 .4],'LineWidth',1.0);
        ylabel('measure of active set','FontSize',10,'FontName',fnt)
        ylim([-0.05 1.2])
        title('computed active set','FontSize',10,'FontName',fnt,'FontWeight','normal')
        set(legend({'\itP\rm_{1}','\itP\rm_{2}'},'Location','best', ...
            'FontSize',7.5,'FontName',fnt),'ItemTokenSize',[12 6])
    else
        set(gca,'YScale','log')
        ylabel('max artificial gap','FontSize',10,'FontName',fnt)
        title('artificial gap','FontSize',10,'FontName',fnt,'FontWeight','normal')
        set(legend({'\itP\rm_{1}','\itP\rm_{2}'},'Location','best', ...
            'FontSize',7.5,'FontName',fnt),'ItemTokenSize',[12 6])
    end
end
exportgraphics(f,[fig 'gap-sweep.pdf'],'ContentType','vector','BackgroundColor','white'); close(f)
disp('extra figures written');

% ---------- gamma0 sweep and unstructured meshes ----------
if isfile('numerics_gamma_mesh.mat')
T=load('numerics_gamma_mesh.mat');

% gamma0: energy norm error against the Nitsche parameter, s held at h^{2k+1}/4
f=figure('Units','centimeters','Position',[2 2 18 7.5],'Visible','off');
for panel=1:2
    subplot(1,2,panel); d=panel; G=T.GAM.(sprintf('p%d',d));
    hc=gobjects(1,numel(G.lv)); hbad=gobjects(1,1);
    for li=1:numel(G.lv)
        ok=G.OK(li,:) & isfinite(G.H1(li,:));
        hc(li)=loglog(G.g0(ok),G.H1(li,ok),['-' mk{li}],'Color',c(li,:), ...
            'LineWidth',1.3,'MarkerSize',5,'MarkerFaceColor',c(li,:)); hold on
    end
    for li=1:numel(G.lv)
        bad=~G.OK(li,:) & isfinite(G.H1(li,:));
        if any(bad)
            hbad=loglog(G.g0(bad),G.H1(li,bad),'x','Color',[.85 .1 .1], ...
                'LineWidth',1.4,'MarkerSize',7);
        end
    end
    % reference slope gamma0^{1/2}: the regularization error grows with mu = gamma*s
    gg=G.g0(G.g0>=20); li=numel(G.lv);
    ref=G.H1(li,find(G.g0>=20,1))*sqrt(gg/gg(1));
    href=loglog(gg,ref,'k--','LineWidth',1.0);
    grid on, box on, set(gca,'FontSize',9,'FontName',fnt)
    xlabel('\it\gamma\rm_0','FontSize',10,'FontName',fnt)
    if panel==1, ylabel('\itH\rm^{1} error','FontSize',10,'FontName',fnt); end
    lg=arrayfun(@(v) sprintf('\\itn\\rm = %d',v),G.lv,'UniformOutput',false);
    hh=[hc, href]; ll=[lg,{'\it\gamma\rm_0^{1/2}'}];
    if isgraphics(hbad), hh=[hh, hbad]; ll=[ll,{'not converged'}]; end
    set(legend(hh,ll,'Location','best','FontSize',7,'FontName',fnt),'ItemTokenSize',[12 6])
    title(sprintf('\\itP\\rm_{%d}',d),'FontSize',10,'FontName',fnt,'FontWeight','normal')
end
exportgraphics(f,[fig 'gamma-sweep.pdf'],'ContentType','vector','BackgroundColor','white'); close(f)

% unstructured meshes: the mesh itself and the convergence, against the structured runs
D=load('exact_conv_data.mat'); DA=D.DATA;
f=figure('Units','centimeters','Position',[2 2 18 6.2],'Visible','off');
subplot(1,3,1)
[t3,xv,yv]=unstructmesh(16);
triplot(t3,xv,yv,'Color',[.35 .35 .35],'LineWidth',0.3);
axis equal tight, box on, set(gca,'FontSize',8,'FontName',fnt)
xlabel('\itx','FontSize',9,'FontName',fnt); ylabel('\ity','FontSize',9,'FontName',fnt)
title('\itn\rm = 16','FontSize',10,'FontName',fnt,'FontWeight','normal')
for panel=1:2
    subplot(1,3,1+panel); k=panel;
    for ic=1:2
        nm=sprintf('%s%d',char('A'+ic-1),k);
        St=DA.(nm); Un=T.UNS.(nm);
        loglog(St.h,St.H1,'-','Color',c(ic,:),'LineWidth',1.3); hold on
        loglog(Un.h,Un.H1,mk{ic},'Color',c(ic,:),'MarkerSize',6, ...
            'MarkerFaceColor','w','LineWidth',1.3);
    end
    h=DA.(sprintf('A%d',k)).h;
    loglog(h,DA.(sprintf('A%d',k)).H1(end)*(h/h(end)).^k,'k:','LineWidth',1.0);
    grid on, box on, set(gca,'FontSize',9,'FontName',fnt)
    xlabel('\ith','FontSize',10,'FontName',fnt)
    if panel==1, ylabel('\itH\rm^{1} error','FontSize',10,'FontName',fnt); end
    if panel==1      % one legend for both panels; the encoding is the same
        yl=ylim; ylim([yl(1) yl(2)*8]);      % headroom, so the box clears the data
        set(legend({'A, structured','A, unstructured','B, structured', ...
            'B, unstructured',sprintf('\\ith\\rm^{%d}',k)}, ...
            'Location','northwest','FontSize',6.5,'FontName',fnt), ...
            'ItemTokenSize',[12 6])
    end
    title(sprintf('\\itP\\rm_{%d}',k),'FontSize',10,'FontName',fnt,'FontWeight','normal')
end
exportgraphics(f,[fig 'unstructured.pdf'],'ContentType','vector','BackgroundColor','white'); close(f)
disp('gamma and mesh figures written');
end
