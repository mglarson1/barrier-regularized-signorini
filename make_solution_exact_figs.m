% Elevation plots of the P2 solutions for the two exact-solution examples,
% with the mesh shown and the active contact zone marked.
% Output: figures/solution-A.pdf, solution-B.pdf
set_publication_graphics;
outdir=fullfile(fileparts(mfilename('fullpath')),'figures');
if ~isfolder(outdir), mkdir(outdir); end
plot_exact_solution(@exact_A,fullfile(outdir,'solution-A.pdf'));
plot_exact_solution(@exact_B,fullfile(outdir,'solution-B.pdf'));
disp('exact-solution figures written')

function plot_exact_solution(exfun,fname)
n = 16; deg = 2; pow = 5; gamma0 = 20;
[tri3,xv,yv] = rectmesh(2*n,n,[-1 1],[0 1]);
[tri,xnod,ynod] = addmidnodes(tri3,xv,yv);
TR = triangulation(tri3,xv,yv); nb = TR.neighbors; nb(isnan(nb)) = 0;
nb = [nb(:,3),nb(:,1),nb(:,2)];
p = multipliers_bottom(tri3,nb,xnod,ynod,0);
ind = find(abs(xnod+1)<1e-10 | abs(xnod-1)<1e-10 | abs(ynod-1)<1e-10);
u = solvereg_ex(tri,tri3,xnod,ynod,p,deg,pow,gamma0,exfun,ind);

nv = length(xv);
f = figure('Units','centimeters','Position',[2 2 14 10],'Visible','off');
trisurf(tri3,xv,yv,u(1:nv),'FaceColor','w','EdgeColor',[.45 .45 .45], ...
    'LineWidth',0.25);
hold on
enod = [1,2;2,3;3,1];
zlift = 6e-3*(max(u)-min(u));
for im = 1:size(p,2)
    iel = p(im).iel; inei = p(im).inei;
    iv = tri(iel,:); xc = xnod(tri3(iel,:)); yc = ynod(tri3(iel,:));
    nod3 = tri3(iel,enod(inei,:));
    midnod = tri(iel,3+inei);
    xcut = xnod(nod3); ycut = ynod(nod3);
    dl = sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
    gamma = gamma0/dl;
    nvec = [ycut(2)-ycut(1),xcut(1)-xcut(2)]/dl;
    [fi,fix,fiy,~] = basisfun(deg,mean(xcut),mean(ycut),xc,yc);
    fin = nvec(1)*fix + nvec(2)*fiy;
    w = dot(fi,u(iv)) - dot(fin,u(iv))/gamma;
    ids = [nod3(1), midnod, nod3(2)];
    if w > 0
        plot3(xnod(ids),ynod(ids),u(ids)+zlift,'r-','LineWidth',3);
    else
        plot3(xnod(ids),ynod(ids),u(ids)+zlift,'-','Color',[.2 .2 .2],'LineWidth',1.5);
    end
end
view(-25,30)
fnt = 'Times New Roman';
set(gca,'FontSize',9,'FontName',fnt)
xlabel('\itx','FontName',fnt), ylabel('\ity','FontName',fnt), zlabel('\itu_{h}','FontName',fnt)
axis tight, box on, grid on
exportgraphics(f,fname,'ContentType','image','Resolution',300,'BackgroundColor','white');
close(f)
end
