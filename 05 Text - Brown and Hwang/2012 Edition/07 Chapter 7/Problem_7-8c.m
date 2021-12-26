%This file is for Problem 7.8, Part (c)

function [hpt] = intpolate(x,h,xpt);
 
i=1;
while x(i) < xpt
    i=i+1;
end;
xmin=x(i-1);
xmax=x(i);
if xmax == xmin
    xmax=xmin+1;
end;
if xmax > 25000
    xmax=25000;
end;
if xmin > 25000
    xmin=24999;
end;
if xmax < 1
    xmax=2;
end;
if xmin < 1
    xmin=1;
end;

hpt=(h(xmin)*(xmax-xpt)+h(xmax)*(xpt-xmin))/(xmax-xmin);
