%function gamma=normalization_factor(ell,M,BC_type,radius,period_truncature)
% Returns a factor used to normalize the correlation operator
% - D: Daley length scale
% - M: 'smoothness' parameter of Matern function
% - BC_type can be equal to: 1 -> Dirichlet 
%                            2 -> Neumann
%                            3 -> Periodic
% - radius: radius of the circle
% - period_truncature: number of terms used in the truncated sums
% - ell: matern length parameter
% - gammma: square root of the normalization factor

function Gamma_matrix=normalization_factor(ell,M,BC_type,grid,period_truncature)

%%%%%%%%%%%%%%%%%%%%%%% BUILDING BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Normalization factor on an infinite line 
gamma_line=2*gamma(M)*sqrt(pi)*ell/gamma(M-0.5);

% Values of beta_j
beta_list=ones(1,M);
for j=1:M-1
    beta_list(j+1)=beta_list(j)*2*(M-j)/(j*(2*M-j-1));
end

% Length of the grid
L=grid(end)-grid(1);

% Distance to a boundary of each point
left_distance=abs(grid(1)-grid);
right_distance=abs(grid(end)-grid);
%distance=2*min([left_distance;right_distance]);
distance=2*(min([left_distance;right_distance])+(grid(2)-grid(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% BOUNDARY CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(BC_type)
    
    
    case(1)%%% Dirichlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        mirror_function=matern(distance,ell,M);
    
        correction_factor=1-mirror_function;
        
    case(2)%%% Neumann %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mirror_function=matern(distance,ell,M);
    
        correction_factor=1+mirror_function;
        
       
    case(3)%%%% Periodic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        correction_factor=1; %initialization
        %K is the number of periods taken into account (this is a truncated
        %sum).
        
        for k=1:period_truncature
            j_sum=0;
            for j=0:M-1
                j_sum=j_sum+beta_list(j+1)*((L*k/ell)^j)*...
                    exp(-L*k/ell);
            end
            correction_factor=correction_factor+2*j_sum;
        end 
        
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_corrected=sqrt(gamma_line./correction_factor);

Gamma_matrix=diag(gamma_corrected)  ;

end