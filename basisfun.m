function [fi,fix,fiy,area]=basisfun(deg,x,y,xc,yc)
% Dispatch to P1 or P2 basis functions; xc,yc are the vertex coordinates.
if deg==1
    [fi,fix,fiy,area]=basescalar(x,y,xc,yc);
else
    [fi,fix,fiy,area]=basescalar2(x,y,xc,yc);
end
end
