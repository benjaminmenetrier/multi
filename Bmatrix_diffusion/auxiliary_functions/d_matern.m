function dc=d_matern(r,ell,M)
%function dc=d_matern(r,ell,M)
%This function returns the derivative of a Matern function.
% Inputs:
% - r: distance used as argument of the function. Can be an array.
% - ell: Matern length parameter
% - M: Matern smoothness parameter
% Outputs
% -dc: values of the function, same size as r.



beta=1; %betais intialized as beta1=1
if (M==1)
    dc_sum=-1;
else
    dc_sum=-r/ell;% initialization of the sum
end

for j=2:M-1
    
    beta=beta*2*(M-j)/(j*(2*M-j-1));%beta is updated
    j_term=(j-r/ell).*(r/ell).^(j-1);
    dc_sum=dc_sum+beta*j_term;
end

dc=(1/ell)*exp(-r/ell).*dc_sum;

end