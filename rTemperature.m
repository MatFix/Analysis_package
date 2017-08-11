R = 160e-9;
h = 20e-9;
Ms = 490e3;
Aex = 1.05e-11;
mu0 = 4*pi*1e-7;

kB = 1.38064852e-23;

t = [0 50 80 100 150 200];

gamma = 1.760859e11;
d = 2*R;
alpha = 0.03;
CD = 2;

fRes = 337e6;
omega = fRes*2*pi;

%%
%pot = 2*pi*mu0*Ms^2*h^2/R*(log(8*R/h) - 1/2);
pot = 9.98*pi*mu0*Ms^2*h^2/R;
r = @(T) sqrt(2*kB*T/pot);

figure
scatter(t,r(t))
hold on
%%
G0 = -2*pi*Ms*mu0*d/(gamma);
lexc = sqrt(2*Aex/(mu0*Ms^2));
rcore = 0.68*lexc*(h/lexc)^(1/3);
%D = -abs(G0)/2*log(R/(rcore*exp(-CD)));
D = -mu0/gamma*Ms*pi*h*log(R/rcore);

r = @(T) sqrt(2*kB*T*abs(G0)/(omega*(G0^2 + alpha^2*D^2)));

%%

scatter(t,r(t))