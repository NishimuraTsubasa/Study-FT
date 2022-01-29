% Black's price under log-normal model 
function [px] = callbl(F, K, r, sigma,time)

d1    = (log(F/K) + ((sigma^2)/2)*time)/(sigma*sqrt(time));
d2    = d1 - sigma*sqrt(time);

% Price of Call option  
px    = exp(-r*time)*(F*normcdf(d1) - K*normcdf(d2)); 

return; 

end

function ncdf = normcdf(x)

a1    =  0.319381530;
a2    = -0.356563782; 
a3    =  1.781477937;
a4    = -1.821255978;
a5    =  1.330274429;
gamma = 0.2316419;

k      = 1.0/(1 + gamma*x); 

if (x >= 0) 
  ncdf   = 1 - normal(x)*(a1*k + a2*k^2 + a3*k^3 + a4*k^4 + a5*k^5); 
else
  ncdf   = 1 - normcdf(-x);  
end

return;

end

function nrml = normal(x) 

nrml = (1.0/sqrt(2.0*pi))*exp(-((x^2)/2.0));

return; 

end