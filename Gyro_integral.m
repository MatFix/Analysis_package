clearvars

delta = 2;
alpha = 1;

k = delta/alpha; 
rho = 0:0.01:3;
beta = 0:0.001:pi;
[R,B] = meshgrid(rho,beta);

costheta = sqrt((alpha^2 * cos(B).^2 + delta^2 * sin(B).^2).^2 + R.^4)./(1 + R.^2);

[GRADCOS,~] = gradient(-costheta,rho,beta);

% f = @(alpha,delta) 4 *(-R.^3 ./sqrt((alpha^2 * cos(B).^2 + delta^2 * sin(B).^2).^2 + R.^4) ./(1 + R.^2) + ...
%                        sqrt((alpha^2 * cos(B).^2 + delta^2 * sin(B).^2).^2 + R.^4)./(1 + R.^2).^2 .* R); 

PHI = @(k) atan2(-k^2 *cos(B),sin(B)) + pi/2;
[~,GRADPHI] = gradient(PHI(k),rho,beta);

% QTY = f(alpha,delta).*GRADPHI;
QTY = 2 * GRADCOS.*GRADPHI;
INTEGRAL = trapz(trapz(QTY))/(length(beta)*length(rho));
disp(INTEGRAL)