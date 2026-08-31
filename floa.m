function f=floa(xm,ym)
%f=(xm-0.5)^2*(xm-0.6);
f = -1000*sin(2*pi*xm);
if(ym < 1/2) f=-1;
else
    f=1;
end
end
