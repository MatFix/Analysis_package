mu0 = 4*pi*1e-7;
Ms = 490e3;
alpha = 0.03;
kB = 1.380648e-23;
t = 20e-9;
Rmax = 160e-9;
Aex = 1.05e-11;
lexc = sqrt(2*Aex/(mu0*Ms^2));
rcore = 0.68*lexc*(t/lexc)^(1/3);
gammaLL = 1.761e11;
G0 = 2*pi*Ms*mu0*t/gammaLL;
deltaT = 5; % K

tau = deltaT/(2*Rmax);
T0 = deltaT/2;

DG = -2*abs(G0); % approximation

mOmSq = 3*G0*1e9; % approximation

FrecX = @(x,y) -mOmSq*x;
FrecY = @(x,y) -mOmSq*y ;
FthX = -(rcore/Rmax)^2*kB*tau;
FthY = -(rcore/Rmax)^2*kB*tau;


% FthX = @(x) - sqrt(2/(gammaLL*dt)*alpha*Ms*kB*t*pi*...
%     ((2 + log(Rmax/rcore))*T0 - (9/4 + log(Rmax/rcore))*tau*x));
% FthY = @(x) - sqrt(2/(gammaLL*dt)*alpha*Ms*kB*t*pi*...
%     ((2 + log(Rmax/rcore))*T0  - (11/4 + log(Rmax/rcore))*tau*x));

% FthX = @(x) 0;
% FthY = @(x) 0;

x = 1e-10;
y = -1e-10;

const = 1/(G0^2 + DG^2 * alpha^2);

dt = 500e-13;
Tstop = 20e-9;

Nstep = Tstop/dt;

Xv = zeros(Nstep+1,1);
Yv = zeros(Nstep+1,1);
time = 0:dt:Tstop;

Xv(1) = x;
Yv(1) = y;

for ii = 1:Nstep
  a = randn(1);
  b = randn(1);
  x = x - const*dt*((FrecX(x,y)+FthY)*alpha*DG + G0*(FrecY(x,y) +FthY)); 
  y = y - const*dt*((FrecY(x,y)+FthY)*alpha*DG - G0*(FrecX(x,y) +FthY)); 
  
  Xv(ii+1) = x;
  Yv(ii+1) = y;
end

%%
figure
plot(time,Xv)
hold on
plot(time,Yv)
legend('x','y')
legend('boxoff')
figure
plot(Xv,Yv)
axis equal
   
    