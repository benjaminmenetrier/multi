%function D=matern_to_daley(ell,M,grid,BC_type,period_truncature)
% This function converts a Matern length parmater into a Daley length scale
% for the diffusion in 1D with different boundary conditions.
% BC_type= 1: Dirichlet
%          2: Neumann
%          3: Periodic


function D=matern_to_daley(ell,M,grid,BC_type,period_truncature)


L=grid(end)-grid(1); %length of the interval

% Distance to a boundary of each point
left_distance=abs(grid(1)-grid);
right_distance=abs(grid(end)-grid);
%distance=2*min([left_distance;right_distance]);
distance=2*(min([left_distance;right_distance])+(grid(2)-grid(1)));


switch(BC_type)
    case(1)% Dirichlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numerator=matern(0,ell,M)-matern(distance,ell,M);
        denominator=d2_matern(0,ell,M)-d2_matern(distance,ell,M);
        
    case(2)% Neumann %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        numerator=matern(0,ell,M)+matern(distance,ell,M);
        denominator=d2_matern(0,ell,M)+d2_matern(distance,ell,M);
        
    case(3) % Periodic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        numerator=matern(0,ell,M); %Initialization of the numerator
        denominator=d2_matern(0,ell,M);%Initialization of denominator
        
        
        %The inifinite sums are truncated. Usually, only the two or the
        %  three firstterms are significant
        for k=1:period_truncature
            numerator=numerator+2*matern(k*L,ell,M); %
            denominator=denominator+2*d2_matern(k*L,ell,M);
        end
 
end



D=sqrt(-numerator/denominator);% Daley length scale


end