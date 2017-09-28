function [ delta_f_FWHM ] = fwhm( X,f )
%% Full Width Half Maximum calculation
half_max = max(X)/2;

A = find(X>half_max);

x1 = f(A(1)-1);
y1 = X(A(1)-1);
x2 = f(A(1));
y2 = X(A(1));
x3 = f(A(end));
y3 = X(A(end));
x4 = f(A(end)+1);
y4 = X(A(end)+1);

delta_f_FWHM = abs(findf(x1,y1,x2,y2,half_max) - findf(x3,y3,x4,y4,half_max));

end

function f = findf(x1,y1,x2,y2,halfmax)
%% linear fit to find intersection (better precision)    
f = halfmax*(x1 - x2)/(y1 - y2) - (x1*y2 - x2*y1)/(y1 - y2);
end