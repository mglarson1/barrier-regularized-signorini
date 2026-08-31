% Elevation plot of the P2 solution with the mesh shown and the active
% contact zone marked. Output: figures/solution-p2.pdf
set_publication_graphics;
n = 32; pow = 5;              % s = h^5/4, the threshold choice for P2
gamma0 = 20;
[tri6,tri3,xnod,ynod] = unitsquaremesh2(n);
TR = triangulation(tri3,xnod(1:max(tri3(:))),ynod(1:max(tri3(:))));
nb = TR.neighbors; nb(isnan(nb)) = 0; nb = [nb(:,3),nb(:,1),nb(:,2)];
p = multipliers(tri3,nb,xnod,ynod);
u = solvereg_eps2(tri6,tri3,xnod,ynod,p,pow,gamma0);

f = figure('Units','centimeters','Position',[2 2 14 10],'Visible','off');
% white faces, mesh edges visible; no elevation coloring
trisurf(tri3,xnod(1:max(tri3(:))),ynod(1:max(tri3(:))),u(1:max(tri3(:))), ...
    'FaceColor','w','EdgeColor',[.45 .45 .45],'LineWidth',0.25);
hold on

% classify each contact edge: active if P(u) = u - sigma_n(u)/gamma > 0
% at the edge midpoint
enod = [1,2;2,3;3,1];
zlift = 4e-3*(max(u)-min(u));  % lift lines slightly above the surface
for im = 1:size(p,2)
    iel = p(im).iel; inei = p(im).inei;
    iv = tri6(iel,:); xc = xnod(tri3(iel,:)); yc = ynod(tri3(iel,:));
    nod3 = tri3(iel,enod(inei,:));
    midnod = tri6(iel,3+inei);
    xcut = xnod(nod3); ycut = ynod(nod3);
    dl = sqrt((xcut(1)-xcut(2))^2+(ycut(1)-ycut(2))^2);
    gamma = gamma0/dl;
    nvec = [ycut(2)-ycut(1),xcut(1)-xcut(2)]/dl;
    xm = mean(xcut); ym = mean(ycut);
    [fi,fix,fiy,~] = basescalar2(xm,ym,xc,yc);
    fin = nvec(1)*fix + nvec(2)*fiy;
    w = dot(fi,u(iv)) - dot(fin,u(iv))/gamma;
    ids = [nod3(1), midnod, nod3(2)];
    if w > 0
        plot3(xnod(ids),ynod(ids),u(ids)+zlift,'r-','LineWidth',3);
    else
        plot3(xnod(ids),ynod(ids),u(ids)+zlift,'-','Color',[.2 .2 .2],'LineWidth',1.5);
    end
end
view(55,28)
set(gca,'FontSize',9,'FontName','Times New Roman')
fnt = 'Times New Roman';
xlabel('\itx','FontName',fnt), ylabel('\ity','FontName',fnt), zlabel('\itu_{h}','FontName',fnt)
axis tight, box on, grid on
outdir=fullfile(fileparts(mfilename('fullpath')),'figures');
if ~isfolder(outdir), mkdir(outdir); end
exportgraphics(f, fullfile(outdir,'solution-p2.pdf'), 'ContentType','image', ...
    'Resolution',300,'BackgroundColor','white');
close(f)
disp('solution figure written')
