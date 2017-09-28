function [K] = Kcalc( eps )
%K calculator for Ni
if norm(eps) > 1
    eps = eps*1e-6;
    disp('eps converted to e-6 form')
end

l100 = -46e-6;
l111 = -24e-6;
c44  = 1.18e11;

ls = 2/5*l100 + 3/5*l111;

K =  - 3/2 * ls * c44*eps;
end

