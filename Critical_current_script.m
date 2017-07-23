% Ni + natural

Msat = 480e3; % Gauss
Aex = 10.5e-12;
alpha = 0.03;

mu0 = 4*pi*1e-7;
hbar = 1.054*1e-34;
qel = 1.602e-19;

gamma = 1.761e11;           %gyromagnetic ratio

% Specific dot
Omega = 220e6*2*pi;

R = 400e-9; % radius
h = 20e-9;  % thickness

pol = 0.4;  % degree of polarization

%%

omega0 = gamma*mu0*Msat;
Omega_tilde = Omega/omega0;

lex = sqrt(2*Aex/(mu0 * Msat^2));

nu = pi*alpha*log(R/lex);
A = (4*pol^(3/2)) /(3*(1 + pol)^3 - 16 * pol^(3/2));

%%

Jcr = 2*nu*Omega_tilde/A * (mu0 * Msat^2 * qel * h / hbar);

Jsw = 1.45*Jcr;

str = sprintf('\nSwitching current density: %.2d A/m^2\n', Jsw);

disp(str)
