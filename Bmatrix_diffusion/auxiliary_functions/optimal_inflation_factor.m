% function mu=optimal_inflation(param_B,param_R,param_H,H_handles)
% Returns the optimal inflation parameter for the variances of a diagonal R
% In order to speed up the process, and because an extreme precision is not
% reuqired, the optimal value is chosen from a small set of possible values

function mu=optimal_inflation_factor(param_B,param_R,param_H,...
    H_handles,resolution,xb,y,xt)


%------------------------- PARAMETERS SETTING ----------------------------%
mu=1;
param_R('inflation_factor')=mu;
param_R('D')=0;
%-------------------------- SOLVERS OPTIONS ------------------------------%
param_solver=containers.Map;
param_solver('precision')=10^(-4);
param_solver('max_iteration_number')=1.2*length(param_B('grid'));
param_solver('solver_type')=2; %1 for basic CG ; 2 for preconditioned CG ;
% 3 for backslash
param_solver('re_orthogonalization')=1; % if equal to 1, the residuals will
%be reorthogonalized at each iteration to compensate the loss of
%orthogonality due to finite precision.
%-------------------------------------------------------------------------%


%-------------------------- DISPLAY OPTIONS ------------------------------%
param_display=containers.Map;
param_display('fontsize')=20;

param_display('uncorrelated_comparison')=0; %If equal to 1, each time a
% quantity is plotted, the same quantity will be plotted in the case where
% R is diagonal as a comparison
param_display('convergence_plot')=0; % Indicators of convergence rate will 
                        % be saved at each iteration of the CG and plottted
%param_display('spectrum')=1;% display the eigenvalue distribution of the 
% hessian matrix
param_display('A_norm_with_bounds')=0; % display a figure with the A norm
%of the error and a set of bounds chosen belox
                        
%----------------------- Bounds on A norm of the error -------------------%
 param_display('kappa_bound')=0;% bound on the A_norm based on the condition
% % number of S.

param_display('spectrum_plot')=0;
%-------------------------------------------------------------------------%




%----------------------- Analysis error computation ----------------------%
    [xa,~,~]=linear_3DVar(param_B,param_R,param_H,H_handles,...
    param_solver,param_display,xb,y,xt);
    old_analysis_error=norm(xa-xt,2);
%-------------------------------------------------------------------------%
found=0;

while(found==0 & mu<15)
    mu=mu+resolution;
    param_R('inflation_factor')=mu;
    %--------------------- Analysis error computation --------------------%
    [xa,~,~]=linear_3DVar(param_B,param_R,param_H,H_handles,...
    param_solver,param_display,xb,y,xt);
    new_analysis_error=norm(xa-xt,2);
    %---------------------------------------------------------------------%
    if (new_analysis_error>old_analysis_error)
        mu=mu-resolution;
        found=1;
    else
        old_analysis_error=new_analysis_error;
    end

end

if(mu>=15)
    disp('stopped optimal_inflation_factor because mu=15')
    disp("Are you sure there isn't a problem ?")
end
param_R('inflation_factor')='optimal';
end