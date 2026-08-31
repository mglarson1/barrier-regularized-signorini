% Figures replacing the tables of Section 6 in the scalar paper.
% The Newton figure is built by make_extra_figs.m, which also holds the
% residual histories it was merged with.
set_publication_graphics;
fnt='Times New Roman';
fig=[fullfile(fileparts(mfilename('fullpath')),'figures') filesep];
if ~isfolder(fig), mkdir(fig); end
c=lines(6);

% ---------- convergence against the exact solutions ----------
D=load('exact_conv_data.mat'); DA=D.DATA;
f=figure('Units','centimeters','Position',[2 2 18 7.5],'Visible','off');
for panel=1:2
    subplot(1,2,panel); deg=panel;
    A=DA.(sprintf('A%d',deg)); B=DA.(sprintf('B%d',deg));
    h=A.h;
    hl(1)=loglog(h,A.H1,'-o','Color',c(1,:),'LineWidth',1.3,'MarkerSize',5,'MarkerFaceColor',c(1,:)); hold on
    hl(2)=loglog(h,A.L2,'--o','Color',c(1,:),'LineWidth',1.0,'MarkerSize',4);
    hl(3)=loglog(B.h,B.H1,'-s','Color',c(2,:),'LineWidth',1.3,'MarkerSize',5,'MarkerFaceColor',c(2,:));
    hl(4)=loglog(B.h,B.L2,'--s','Color',c(2,:),'LineWidth',1.0,'MarkerSize',4);
    hr=loglog(h,A.H1(end)*(h/h(end)).^deg,'k-','LineWidth',1.0);
    if deg==2
        loglog(B.h,B.H1(end)*(B.h/B.h(end)).^1.5,'k:','LineWidth',1.2);
    end
    grid on, box on, set(gca,'FontSize',9,'FontName',fnt)
    xlabel('\ith','FontSize',10,'FontName',fnt)
    if panel==1, ylabel('error','FontSize',10,'FontName',fnt); end
    lg={'A, \itH\rm^{1}','A, \itL\rm^{2}','B, \itH\rm^{1}','B, \itL\rm^{2}', ...
        sprintf('\\ith\\rm^{%d}',deg)};
    if deg==2, lg=[lg,{'\ith\rm^{3/2}'}]; end
    yl=ylim; ylim([yl(1) yl(2)*10]);     % headroom for the legend
    set(legend(lg,'Location','northwest','FontSize',7.5,'FontName',fnt),'ItemTokenSize',[12 6])
    title(sprintf('\\itP\\rm_{%d},  \\its\\rm = \\ith\\rm^{%d}/4',deg,2*deg+1),'FontSize',10,'FontName',fnt,'FontWeight','normal')
end
exportgraphics(f,[fig 'exact-convergence.pdf'],'ContentType','vector','BackgroundColor','white'); close(f)

% ---------- effect of the regularization parameter ----------
SV=D.SV; pw=[3 5 7];
f=figure('Units','centimeters','Position',[2 2 18 7.5],'Visible','off');
for panel=1:2
    subplot(1,2,panel)
    for k=1:3
        Q=SV.(sprintf('p%d',pw(k)));
        y=Q.H1; if panel==2, y=Q.L2; end
        loglog(Q.h,y,'-o','Color',c(k,:),'LineWidth',1.3,'MarkerSize',5, ...
            'MarkerFaceColor',c(k,:)); hold on
    end
    Q=SV.p7;
    if panel==1, loglog(Q.h,Q.H1(end)*(Q.h/Q.h(end)),'k-','LineWidth',1.0);
    else, loglog(Q.h,Q.L2(end)*(Q.h/Q.h(end)).^2,'k-','LineWidth',1.0); end
    grid on, box on, set(gca,'FontSize',9,'FontName',fnt)
    xlabel('\ith','FontSize',10,'FontName',fnt)
    lg=[arrayfun(@(a) sprintf('\\its\\rm = \\ith\\rm^{%d}/4',a),pw,'UniformOutput',false)];
    if panel==1
        ylabel('\itH\rm^{1} error','FontSize',10,'FontName',fnt); lg=[lg,{'\ith'}];
        title('energy norm','FontSize',10,'FontName',fnt,'FontWeight','normal')
    else
        ylabel('\itL\rm^{2} error','FontSize',10,'FontName',fnt); lg=[lg,{'\ith\rm^{2}'}];
        title('\itL\rm^{2} norm','FontSize',10,'FontName',fnt,'FontWeight','normal')
    end
    yl=ylim; ylim([yl(1) yl(2)*10]);     % headroom for the legend
    set(legend(lg,'Location','northwest','FontSize',7.5,'FontName',fnt),'ItemTokenSize',[12 6])
end
exportgraphics(f,[fig 'sparam-effect.pdf'],'ContentType','vector','BackgroundColor','white'); close(f)
disp('scalar figures written')
