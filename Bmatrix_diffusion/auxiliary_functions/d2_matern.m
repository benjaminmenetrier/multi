function d2c=d2_matern(r,ell,M)
%function d2c=d2_matern(r,ell,M)
%This function returns the second order derivative of a Matern function.
% Inputs:
% - r: distance used as argument of the function. Can be an array.
% - ell: Matern length parameter
% - M: Matern smoothness parameter
% Outputs
% -d2c: values of the function, same size as r.



beta=1; %betais intialized as beta1=1
if (M==1)
    d2c_sum=-1;
else
    d2c_sum=1-r/ell;% initialization of the sum
end

for j=2:M-1
    
    beta=beta*2*(M-j)/(j*(2*M-j-1));%beta is updated
    j_term=( (r/ell).^(j-2) ) .* ( j*(1-j) +2*j*(r/ell) -(r/ell).^2 ); 
    d2c_sum=d2c_sum+beta*j_term;
end

d2c=-(1/ell^2)*exp(-r/ell).*d2c_sum;

end