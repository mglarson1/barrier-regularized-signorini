% Generate convergence figures for Section 6 from the sweep results.
% Output: figures/regsweep-p1.pdf, figures/regsweep-p2.pdf
set_publication_graphics;
outdir = [fullfile(fileparts(mfilename('fullpath')),'figures') filesep];
if ~isfolder(outdir), mkdir(outdir); end
S = load('results_p1.mat');
plot_sweep(S.levels, S.pows, S.H1, [outdir 'regsweep-p1.pdf']);
S = load('results_p2.mat');
plot_sweep(S.levels, S.pows, S.H1, [outdir 'regsweep-p2.pdf']);
disp('figures written')

function plot_sweep(levels, pows, H1, fname)
h = 1./levels(:);
f = figure('Units','centimeters','Position',[2 2 12 9],'Visible','off');
fnt = 'Times New Roman';
mk = {'o','s','^','d','v'};
cols = lines(numel(pows));
hl = gobjects(1,numel(pows));
for pi = 1:numel(pows)
    hl(pi) = loglog(h, H1(:,pi), ['-' mk{pi}], 'Color', cols(pi,:), ...
        'LineWidth', 1.1, 'MarkerFaceColor', cols(pi,:), 'MarkerSize', 5);
    hold on
end
% dashed reference lines with the predicted slope (alpha-1)/2,
% anchored at the finest-mesh data point
for pi = 1:numel(pows)
    r = (pows(pi)-1)/2;
    loglog(h, H1(end,pi)*(h/h(end)).^r, '--', 'Color', [.6 .6 .6], 'LineWidth', .8);
end
grid on, box on
set(gca,'FontSize',10,'FontName',fnt)
xlabel('\ith','FontSize',11,'FontName',fnt)
ylabel('|\itu_{h,s}\rm - \itu_{h,0}\rm|_{\itH\rm^{1}(\Omega)}','FontSize',11,'FontName',fnt)
leg = arrayfun(@(a) sprintf('\\its\\rm = \\ith\\rm^{%d}/4', a), pows, 'UniformOutput', false);
yl=ylim; ylim([yl(1) yl(2)*100]);        % headroom for the legend
set(legend(hl, leg{:}, 'Location','northwest', 'FontSize', 9, 'FontName', fnt),'ItemTokenSize',[12 6])
xlim([min(h)*0.8, max(h)*1.25])
exportgraphics(f, fname, 'ContentType', 'vector', 'BackgroundColor', 'white');
close(f)
end
