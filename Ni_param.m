% Nickel parameters
% (Gilbert 2016)

E = 2e11;

c11 = 2.5e11;
c12 = 1.6e11;
c44 = 1.18e11;

l100 = -46e-6;
l111 = -24e-6;

Aex = 1.05e-11;
alpha = 0.05;
Ms = 4.8e5;


c = c12*ones(3,3) - diag(c12*ones(1,3)) + diag(c11*ones(1,3));

sx = 0.1e9;
sy = 0;
% 
% stress = [sx 0 0; 0 sy 0 ; 0 0 0];
% 
% eps = c\stress;

ls = 2/5*l100 + 3/5*l111;

eps = @(stress) c11\stress;

str = (0.2:0.2:5)*1e9;

B1 = -3 * ls *c44;

epsxx = linspace(0,1500,50)*1e-6;
epsyy = 0;

Ku1 = - 3/2 * ls * c44 * (epsxx - epsyy);
 
plot(Ku1)