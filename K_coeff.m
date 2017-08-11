function [] = K_coeff(R)

R = R/2;

% R = 160e-9;
h = 20e-9;
Ms = 490e3;
Aex = 1.05e-11;
mu0 = 4*pi*1e-7;

gamma = 1.760859e11;
d = 2*R;
alpha = 0.03;
CD = 2;

fRes = 337e6;
omega = fRes*2*pi;

%%

lexc = sqrt(2*Aex/(mu0*Ms^2));
rcore = 0.68*lexc*(h/lexc)^(1/3);

Acoeff = Aex*log(R/rcore);
Kcoeff = @(Ku1) Ku1*(R^2);

c = @(Ku1) sqrt(1./(1 - Kcoeff(Ku1)./(4*Acoeff)));

%%

Ku1 = 0:0.1:15;

%%

RR = R^2/2;
BB = @(Ku1) RR*Ku1;

cc = @(Ku1) 1 + 1/3 - BB(Ku1)/(24*Acoeff) - sqrt(64*Aex^2 + 32*Aex*BB(Ku1) + BB(Ku1).^2)./(24*Aex);

plot(Ku1,cc(Ku1*1e3))

end