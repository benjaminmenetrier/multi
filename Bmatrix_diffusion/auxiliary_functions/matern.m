function c=matern(r,ell,M)
%function c=matern(r,ell,M)
%This function returns a Matern function.
% Inputs:
% - r: distance used as argument of the function. Can be an array.
% - ell: Matern length parameter
% - M: Matern smoothness parameter
% Outputs
% -c: values of the function, same size as r.


beta=1; %beta is initialized as beta0=1;
c_sum=1;% initialization of the sum in the expression of c


    for j=1:M-1
        beta=beta*2*(M-j)/(j*(2*M-j-1));%beta is updated
        c_sum=c_sum+beta*((r/ell).^j);
    end
c=c_sum.*exp(-r/ell);
    
end