function [r, J] = vasicek_ex1r(x,par)
%%%%%%%%%%%%%%%%%%%%%   
 
load vasicek_data 
discountFactors = vasicek; 
r0              = shortRate;

theta           = x(1); 
alpha           = x(2); 
sigma           = x(3); 

for i = 1:length(discountFactors),
    A(i) = (1.0/alpha) * (1 - exp(-alpha*i)); 
    V(i) = (sigma^2/(2*alpha^3))*(-3 - exp(-2*alpha*i) + 4*exp(-alpha*i) + 2*alpha*i);          
    B(i) = exp(-r0*A(i) - (theta/alpha)*(i - A(i)) + 0.5*V(i));  
end

% Residuals 
r = (discountFactors - B)'; 

end