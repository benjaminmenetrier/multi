addpath('operators');
addpath('auxiliary_functions');
%----------------------------- GEOMETRY ----------------------------------%
N=1000; % Number of points on the domain

L=1000; %length of the domain

BC_type=3;% 1:Dirichlet, 2:Neumann, 3:Periodic

ratio=2;

hb=L/(N-1);

model_grid=[-L/2:hb:L/2]';
%-------------------------------------------------------------------------%


%------------------------- DIFFUSION OPERATORS ---------------------------%
param_B=containers.Map;

param_B('grid')=model_grid;

param_B('sigma')=0.8;  %standard deviation of the background

param_B('D')=60;  %correlation length of the background

param_B('M')=8;   %roughness parameter of the background

param_B('period_truncature')=10;%Number of periods taken into account.
%Usually 2 or 3 are enough.

param_B('inflation_factor')=1; %for no variance inflation, set it to 1

param_B('BC_type')=BC_type; 
%-------------------------------------------------------------------------%

B=diffusion_operator(param_B,eye(N));
B_inv=diffusion_inverse(param_B,eye(N));
B_sqrt=diffusion_sqrt(param_B,eye(N));
