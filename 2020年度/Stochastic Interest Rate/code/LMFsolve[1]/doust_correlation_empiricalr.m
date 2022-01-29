function r = doust_correlation_empiricalr(x,rho)
%%%%%%%%%%%%%%%%%%%%%   

%load rho
%load rho_Mar2011
load rho_Mar2013

N          = 20;       % number of forward rates 
delta_t    = 1;        % accural period 

% two parameter Doust:
Beta   = x(1) * (1./((1:(N-1)).^x(2)));

% five parameter Doust
%Beta   = x(1) + x(2) * (1./((1:(N-1)).^1)) + x(3) * (1./((1:(N-1)).^2)) ... 
%       + x(4) * (1./((1:(N-1)).^3)) + x(5) * (1./((1:(N-1)).^4));
a      = exp(-Beta*delta_t);

sum = 1;
for i = 1:(length(Beta)+1),      
  for j = 1:(length(Beta)+1),
     
     if (j < i), 
       for k = j:(i-1),
         sum = sum * a(k);
       end 
       g(i,j) = sum; 
    
       g(j,i) = g(i,j);
       sum = 1;
     end   
  end
  g(i,i) = 1; 
end

% Residuals
k = 1; 
for i = 1:N, 
   for j = (i+1):N,
     r(k,1) = [g(i,j) - rho(i,j)]; 
     k = k + 1; 
   end
end

end

