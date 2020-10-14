%function ell=daley_to_matern(D,M,BC_type,radius,period_truncature)
% Returns the Matern parameter ell corresponding to D and M, by using a
% dichotomy solver with the function matern_to_daley
% - D: Daley length scale
% - M: 'smoothness' parameter of Matern function
% - BC_type can be equal to: 1 -> Dirichlet
%                            2 -> Neumann
%                            3 -> Periodic
% - L: length of the domain
% - period_truncature: number of terms used in the truncated sums
% - ell: matern length parameter
function ell=daley_to_matern(D,M,grid,BC_type,period_truncature)

% fzero uses dichotomy to find ell
interval=[10^(-16) 10*(grid(end)-grid(1))];
ell=fzero(@fun,interval);


    function out=fun(ell)
        out=matern_to_daley(ell,M,grid,BC_type,period_truncature)-D;
    end
end