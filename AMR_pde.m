%% AMR calculator
%  
%  Given .ovf input files giving the magnetization at each cell, it solves
%  for each file the nonlinear Laplace equation
%  
%  div * (sigma(grad(u)) * grad(u)) = 0
% 
%  where sigma is the AMR conductivity, depending on the reciprocal
%  directions of the reduced magnetization m and the current density at any
%  given point according to
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

global xx yy Mx My rCoeff

%% Simulation parameters

Nx = 40;
Ny = 40;
c = 5e-9;                          % cell size
t = 30e-9;                         % device thickness

% Device resistivity parameters:
%
% rhoTot - average resistivity
% rhoT - transveral resistivity
% rhoL - longitudinal resistivity
% deltaRho/rhoTot - AMR ratio
%
% - data from from Bogart and Atkinson (2009) and McGuire and Potter (1975)
% - valid for Py with Ni80Fe20

t_nm = t*1e9;
rhoTot = 1.2935*t_nm^(-0.4716) * 1e-6;  % Ohm.m
deltaRho = 0.03*rhoTot;
rhoT = (3*rhoTot - deltaRho)/3;

% this coefficient goes into the cfunction calculation
rCoeff = deltaRho/rhoT;

% files renaming
rename = false;

xx = linspace(-1,1,Nx);
yy = linspace(-1,1,Ny);

%% Files folder and files rename

dailyFolder = 'D:\Program Files\mumax\Simulazioni\AMR\';
simulationFolder = 'nanodots_hysteresis_1_dot_longer\';

folder = [dailyFolder simulationFolder];            % folder containing files
PythonScript = 'batchRenamer.py';                   % Python rename script

% calls the MATLAB wrapper to a pyhton function (requires python 3.5 or
% above) which renames files according to the scheme A_ii.ovf, ii = 1, ii++

if rename
    renameFiles(folder,PythonScript);
end

%%

fid = fopen([folder 'table.txt']);

% skip first line
fgets(fid);                                     
matrix = cell2mat(textscan(fid, '%f %f %f %f %f %f %f%*[^\n]'));
fclose(fid);
B = matrix(:,5);
clear matrix

N = length(B);

%% preallocation

R = zeros(size(B));
n = 0;

%% Create a PDE Model with a single dependent variable

model = createpde(1);

% PDE factors
a = 0;
f = 0;

load('variables')

g = decsg(gd,sf,ns);
geometryFromEdges(model,g);

% potential difference for boundary conditions
deltaV = 2;         % Volts

applyBoundaryCondition(model,'Edge',5, 'u', +deltaV/2);
applyBoundaryCondition(model,'Edge',6, 'u', +deltaV/2);
applyBoundaryCondition(model,'Edge',7, 'u', -deltaV/2);
applyBoundaryCondition(model,'Edge',8, 'u', -deltaV/2);

% generate the mesh according the the model described
generateMesh(model,'Hmax',0.2);

% grid for the interpolation of the solution
[X,Y] = meshgrid(xx);
querypoints = [X(:),Y(:)]';

% line of integration in the R determination
lx = 20;

%% AMR calculation

for kk = 1:N
    fid = fopen([folder 'A_' num2str(kk) '.ovf']);
    
    % read the Mx, My columns
    Mread = cell2mat(textscan(fid,'%f %f %*[^\n]','CommentStyle','#'));
    
    fclose(fid);
    
    mx = Mread(1:(Nx*Ny),1);
    my = Mread(1:(Nx*Ny),2);
    
    % Mx and My are arranged so that their indices M(m,n) go like y,x:
    % mx = [1 2 3 4 5 6 7 8]'
    %
    % Mx = [5 6 7 8;
    %       1 2 3 4]
    % This ensures that the PDE results and the matrices are arranged in
    % the same way
    
    Mx = flip((reshape(mx, [Ny,Nx]))');
    My = flip((reshape(my, [Ny,Nx]))');
 
    clear mx my Mread
    
    % solves the non-linear laplace equation using the coefficient given by
    % cfunction(x,y,ux,uy)
    
    u = pdenonlin(model,'cfunction(x,y,ux,uy)',a,f);
    
    % interpolate the results on the cartesian grid
    
    result = createPDEResults(model,u);
    uintrp = interpolateSolution(result,querypoints);
    uintrp = reshape(uintrp,size(X));
    
    % calculate the the gradients of u on the grid
    
    [Fx, Fy] = gradient(uintrp);
    
    Fx(isnan(Fx)) = 0;
    Fy(isnan(Fy)) = 0;
    
    mod_squared = Fx(:,lx).^2 + Fy(:,lx).^2;
    
    if mean(mod_squared) == 0
        mod_squared = 1;
    end
    
    % calculate the effective conductivity at the chosen x coordinate (lx)
    sigma_line = 1./(rhoT + deltaRho./mod_squared.*(Mx(:,lx).*Fx(:,lx) + My(:,lx).*Fy(:,lx)).^2);
    sigma_line(isnan(sigma_line)) = 0;
    
    % the conductance is given by the integral of J and sigma along the y
    % coordinate
    S = - t*trapz(Fx(:,lx).*sigma_line);
    
    % final resistance value
    R(kk) = abs(deltaV/S);
    
    % Print state to console
    fprintf_r('Step completed: %i',kk);
end

fprintf_r('reset');
fprintf('\n');

%% Plotting

% FEM mesh and structure plot
figure
pdemesh(model);
axis equal

title('FEM mesh on the structure')

% direction of variation of B, for the hysteresis plotting
G = gradient(B);
goingUp   = G>=0;
goingDown = G<0;

% AMR plot
 figure
a = plot(B(goingUp),R(goingUp),'b-o','LineWidth',1.5,...
     'MarkerFaceColor','b','MarkerSize',4);
hold on
b = plot(B(goingDown),R(goingDown),'r-o','LineWidth',1.5,...
    'MarkerFaceColor','r','MarkerSize',4);

% legend([a b], 'B ascending', 'B descending')

xlabel('B [T]')
ylabel('R [\Omega]')
title('Anisotropic magnetoresistance with current parallel to applied field')

% electric potential distribution at B = Bmax
figure
pdeplot(model,'xydata',u,'zdata',u)
colormap parula
grid on
string = sprintf('Potential distribution at B = %i mT', max(B)*1e3);
title(string)