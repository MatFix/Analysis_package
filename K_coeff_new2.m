

R = 500e-9;

multiplier = 1;

R = R/multiplier;

% R = 160e-9;
h = 20e-9;
Ms = 490e3;
Aex = 1.05e-11;
mu0 = 4*pi*1e-7;

gamma = 1.760859e11;
alpha = 0.03;
CD = 2;

fRes = 337e6;
omega = fRes*2*pi;

a1 = 1.904;
a2 = -3.406;
a3 = 3.5;

b1 = -0.3082;
b2 = 1.61;
b3 = -0.5507;

p1 = 1.465e-27;
p2 = -1.963e-24;

Ku = -(0:0.1:10)*1e3;

% c = 1:0.001:1.5;

%%

Rmax = R;
A = Aex;

lexc = sqrt(2*Aex/(mu0*Ms^2));
rcore = 0.68*lexc*(h/lexc)^(1/3);

r = linspace(rcore,R,100);

K = -3000;

m = 2*K*r.^2/A;

sol1 = sqrt((-1/2)+(-1/2).*3.^(-1/2).*(3+(-2).*m+m.^2.*(54.*m+m.^3+6.*3.^(1/2) ...
  .*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+(54.*m+m.^3+6.*3.^(1/2).*(27.* ...
  m.^2+m.^4).^(1/2)).^(1/3)).^(1/2)+(-1/2).*(2+(-4/3).*m+(-1/3).* ...
  m.^2.*(54.*m+m.^3+6.*3.^(1/2).*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+( ...
  -1/3).*(54.*m+m.^3+6.*3.^(1/2).*(27.*m.^2+m.^4).^(1/2)).^(1/3)+( ...
  -1/4).*3.^(1/2).*(8+8.*m).*(3+(-2).*m+m.^2.*(54.*m+m.^3+6.*3.^( ...
  1/2).*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+(54.*m+m.^3+6.*3.^(1/2).*( ...
  27.*m.^2+m.^4).^(1/2)).^(1/3)).^(-1/2)).^(1/2));

sol2 = sqrt((-1/2)+(-1/2).*3.^(-1/2).*(3+(-2).*m+m.^2.*(54.*m+m.^3+6.*3.^(1/2) ...
  .*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+(54.*m+m.^3+6.*3.^(1/2).*(27.* ...
  m.^2+m.^4).^(1/2)).^(1/3)).^(1/2)+(1/2).*(2+(-4/3).*m+(-1/3).* ...
  m.^2.*(54.*m+m.^3+6.*3.^(1/2).*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+( ...
  -1/3).*(54.*m+m.^3+6.*3.^(1/2).*(27.*m.^2+m.^4).^(1/2)).^(1/3)+( ...
  -1/4).*3.^(1/2).*(8+8.*m).*(3+(-2).*m+m.^2.*(54.*m+m.^3+6.*3.^( ...
  1/2).*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+(54.*m+m.^3+6.*3.^(1/2).*( ...
  27.*m.^2+m.^4).^(1/2)).^(1/3)).^(-1/2)).^(1/2));

sol3 = sqrt((-1/2)+(1/2).*3.^(-1/2).*(3+(-2).*m+m.^2.*(54.*m+m.^3+6.*3.^(1/2) ...
  .*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+(54.*m+m.^3+6.*3.^(1/2).*(27.* ...
  m.^2+m.^4).^(1/2)).^(1/3)).^(1/2)+(-1/2).*(2+(-4/3).*m+(-1/3).* ...
  m.^2.*(54.*m+m.^3+6.*3.^(1/2).*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+( ...
  -1/3).*(54.*m+m.^3+6.*3.^(1/2).*(27.*m.^2+m.^4).^(1/2)).^(1/3)+( ...
  1/4).*3.^(1/2).*(8+8.*m).*(3+(-2).*m+m.^2.*(54.*m+m.^3+6.*3.^(1/2) ...
  .*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+(54.*m+m.^3+6.*3.^(1/2).*(27.* ...
  m.^2+m.^4).^(1/2)).^(1/3)).^(-1/2)).^(1/2));

sol4 = sqrt((-1/2)+(1/2).*3.^(-1/2).*(3+(-2).*m+m.^2.*(54.*m+m.^3+6.*3.^(1/2) ...
  .*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+(54.*m+m.^3+6.*3.^(1/2).*(27.* ...
  m.^2+m.^4).^(1/2)).^(1/3)).^(1/2)+(1/2).*(2+(-4/3).*m+(-1/3).* ...
  m.^2.*(54.*m+m.^3+6.*3.^(1/2).*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+( ...
  -1/3).*(54.*m+m.^3+6.*3.^(1/2).*(27.*m.^2+m.^4).^(1/2)).^(1/3)+( ...
  1/4).*3.^(1/2).*(8+8.*m).*(3+(-2).*m+m.^2.*(54.*m+m.^3+6.*3.^(1/2) ...
  .*(27.*m.^2+m.^4).^(1/2)).^(-1/3)+(54.*m+m.^3+6.*3.^(1/2).*(27.* ...
  m.^2+m.^4).^(1/2)).^(1/3)).^(-1/2)).^(1/2));



%%
plot(r,abs(sol1))
hold on
plot(r,abs(sol2))
plot(r,abs(sol3))
plot(r,abs(sol4))