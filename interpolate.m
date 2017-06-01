function [ r ] = interpolate( f0,d1,f1,d2,f2 )
% interpolation function

b = (f2 - f1)/(d2 - d1);
a = ((f2-f0)/d2 - (f0-f1)/(-d1)) / (d2 - d1);

r = - b / (2 * a);

end

