

R = 100e-9;

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

K = -(0:0.1:10)*1e3;

% c = 1:0.001:1.5;

%%

lexc = sqrt(2*Aex/(mu0*Ms^2));
rcore = 0.68*lexc*(h/lexc)^(1/3);

A = h*Aex*log(R/rcore)*pi;
B = pi*(R.^2 - rcore.^2)/2*h;

% c = (0.3406E1.*A+(-0.161E1).*B).*(0.3808E1.*A+(-0.6164E0).*B).^(-1);

aa = 1.0075;
bb = 0.4856;

% c = (1 - aa *B/A * K);

c = 0.5E0.*(0.662078E55.*A+(-0.160756E55).*B.*K).^(-1).*(0.125426E56.* ...
  A+(-0.667016E55).*B.*K+(((-0.125426E56).*A+0.667016E55.*B.*K).^2+( ...
  -0.4E1).*(0.662078E55.*A+(-0.160756E55).*B.*K).*(0.592184E55.*A+( ...
  -0.341297E31).*K+(-0.37567E55).*B.*K+0.509424E28.*K.^2)).^(1/2)); ...

% c = 0.5E0.*(0.662078E55.*A+(-0.160756E55).*B.*K).^(-1).*(0.125426E56.* ...
%   A+(-0.667016E55).*B.*K+(-0.1E1).*(((-0.125426E56).*A+0.667016E55.* ...
%   B.*K).^2+(-0.4E1).*(0.662078E55.*A+(-0.160756E55).*B.*K).*( ...
%   0.592184E55.*A+(-0.341297E31).*K+(-0.37567E55).*B.*K+0.509424E28.* ...
%   K.^2)).^(1/2));
plot(abs(K/1e3),c)
hold on