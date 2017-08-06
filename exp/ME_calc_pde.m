%% AMR calculator
%  
%  Given .ovf input files giving the magnetization at each cell, it solves
%  for each file the nonlinear Laplace equation
%  
%  div * (sigma(grad(u)) * grad(u)) = 0
% 
%  where sigma is the AMR conductivity, depending on the relative
%  directions of the reduced magnetization m and the current density J at
%  any given point according to
%
%  1/sigma = rhoT + deltaRho * cos(phi)^2
%
%  where phi is the angle between the two vectors. Once the u has been
%  obtained, the total resistance R can be calculated; the program also
%  displays the AMR chart using the B values extracted from the table.txt
%  file.
%
%  Matteo Fissore, 2017

clear variables
close all

%% Create a PDE Model with a single dependent variable

model = createpde;
load('structure')

g = decsg(gd,sf,ns);
geometryFromEdges(model,g);

stressFactor = 0;
stringStressFactor = num2str(stressFactor,'%10.5e\n');

% PDE factors
a = 0;
c = -1;
% f = [stringStressFactor '*uy./(ux .* sqrt(uy.^2 / ux.^2 + 1))'];
f = [stringStressFactor '*cos(u).*sin(u)'];

boundaryfun = @(region,state) atan2(region.y,region.x) + 2*pi*(region.y<0);
% boundaryfun = @(region,state)region.x.^2;
applyBoundaryCondition(model,'edge',1:8,...
                      'u',boundaryfun,'Vectorized','on');

generateMesh(model,'Hmax',0.02);
u = pdenonlin(model,c,a,f,'U0',1,'Jacobian', 'full','Report', 'on');

%% PLOT
figure
pdeplot(model,'XYData',u,'ZData',u)

figure
xx = -1:0.05:1;
[X,Y] = meshgrid(xx);
querypoints = [X(:),Y(:)]';

result = createPDEResults(model,u);
uintrp = interpolateSolution(result,querypoints);
beta = reshape(uintrp,size(X));

MX = -cos(beta);
MY =  sin(beta);

quiver(X,Y,MX,MY)

