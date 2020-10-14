%function out=diffusion_operator(param_diff,input)
%This function applies a covariance operator based on diffusion,
%parametrized by param_diff.

function output=diffusion_operator(param_diff,input)

grid=param_diff('grid');
sigma=param_diff('sigma');
D=param_diff('D');
M=param_diff('M');
period_truncature=param_diff('period_truncature');
BC_type=param_diff('BC_type');

if (D==0)%If D=0, then the operator just apply variances
    sigma=sigma*param_diff('inflation_factor'); % variances can be inflated
    %to try to compensate the absence of correlations
    output=sigma*input*sigma;
    
else %if D>0, correlations are generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUILDING BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Length of the domain
L=grid(end)-grid(1);

% Matern length scale   
ell=daley_to_matern(D,M,grid,BC_type,period_truncature);
    
% Spatial resolution
h=grid(2)-grid(1);

% Number of points
N=length(grid);

% Normalized Matern length scale
ell_tild=ell/h;

% (I-ell_t^2 Laplacian) discretized with finite differences
T=finite_differences_T(ell_tild,N,BC_type);

% Normalization factor
gamma=normalization_factor(ell,M,BC_type,grid,period_truncature);
%imagesc(gamma)
% figure()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% OPERATORS APPLICATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Operators are applied symmetrically, and can be summed up as:
% output= sigma*gamma*W^(-1/2)*A^(-M)W^(-1/2)*gamma*sigma*input

%variances operator (1/2)
input=sigma*input; 


%normalization operator (1/2)
input =gamma*input;

% Gram inverse matrix (1/2)
input = input/sqrt(h);  


% Finite differences operator
    for k = 1:M
        input = T\input;
    end
% Gram inverse matrix (2/2)
input = input/sqrt(h); 


%normalization operator (2/2)
input =gamma*input;

%variances operator (2/2)
output=sigma*input; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end