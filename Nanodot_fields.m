mu0 = 4*pi*1e-7;

%Ni params
Ms = 480e3;
A = 10.5e-12;

lex = sqrt(2*A/(mu0*Ms^2));

% dot params 
R = 160e-9;
h = 40e-9;

beta = h/R;

%%
t = 0:0.001:500;
K = @(mu,beta) (1 - (1 - exp(-beta*t))./(beta*t))./t .*(besselj(mu,t)).^2;

f1 = K(1,beta);
f2 = K(2,beta);
f1(isnan(f1)) = 0;
f2(isnan(f2)) = 0;

F1 = trapz(t,f1);
F2 = trapz(t,f2);

% Nucleation field [T]

Bn = mu0*Ms*(F1 - F2 - 1/pi *(lex/R)^2);

fprintf('Nucleation field: %f mT\n',Bn*1e3)

% Expulsion field [T]

Be = 2/(4*pi)*mu0*Ms*(2*pi*F1 - 1/2 * (lex/R)^2);

fprintf('Expulsion field: %f mT\n',Be*1e3)