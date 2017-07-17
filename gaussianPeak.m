function [cx,cy,sigmax,sigmay] = gaussianPeak(Mmat,m,n,n_int)
% Gaussian fit
%  gaussianPeak reads the magnetization matrix Mmat and the approximate
%  peak locations m,n, and the number of points n_int (actual number of fit
%  points per dimension n_int*2 + 1); it performs a 2D gaussian fit
%  returning peak coordinates (cx,cy) and standard deviations (sigmax,
%  sigmay).
%
%  For the parent function see also: PEAKFINDER
%

xx = log(abs(Mmat(n,m-n_int:m+n_int)));
yy = log(abs(Mmat(n-n_int:n+n_int,m)))';

Px = polyfit(m-n_int:m+n_int,xx,2);
Py = polyfit(n-n_int:n+n_int,yy,2);

% Px = linFit(m-n_int:m+n_int,xx,2);
% Py = linFit(n-n_int:n+n_int,yy,2);

cx = -Px(2)/(2*Px(1));
cy = -Py(2)/(2*Py(1));

sigmax = sqrt(-1/(2*Px(1)));
sigmay = sqrt(-1/(2*Py(1)));
end