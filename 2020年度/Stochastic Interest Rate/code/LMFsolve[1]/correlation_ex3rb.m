function [r, J] = correlation_ex3rb(x,par)
%%%%%%%%%%%%%%%%%%%%%   

load doust_Mar2013 
%load doust_Mar2011

N = 20;
% partial derivative of f wrt to theta 
d = zeros(N,4,3);

theta = reshape(x,3,N)';
for i = 1:N,
  f(i,:)   = [ cos(theta(i,1)),  cos(theta(i,2))*sin(theta(i,1)),  cos(theta(i,3))*sin(theta(i,1))*sin(theta(i,2)), ... 
                                                                   sin(theta(i,1))*sin(theta(i,2))*sin(theta(i,3))];

  d(i,:,1) = [-sin(theta(i,1)),  cos(theta(i,2))*cos(theta(i,1)),  cos(theta(i,3))*cos(theta(i,1))*sin(theta(i,2)), ... 
                                                                   cos(theta(i,1))*sin(theta(i,2))*sin(theta(i,3))];

  d(i,:,2) = [ 0.0,             -sin(theta(i,2))*sin(theta(i,1)),  cos(theta(i,3))*sin(theta(i,1))*cos(theta(i,2)), ... 
                                                                   sin(theta(i,1))*cos(theta(i,2))*sin(theta(i,3))];
  
  d(i,:,3) = [ 0.0,              0.0,                             -sin(theta(i,3))*sin(theta(i,1))*sin(theta(i,2)), ... 
                                                                   sin(theta(i,1))*sin(theta(i,2))*cos(theta(i,3))];
end
g = f * f';

% Residuals
r = zeros([N*(N-1)/2,1]);
J = zeros([N*(N-1)/2,3*N]);

k = 1;
for i = 1:N, 
   for j = (i+1):N,
     r(k,1) = [g(i,j) - rho(i,j)]; 
     for m = 1:N, 
       for l = 1:3, 
         if (m==j),
           J(k,3*(m-1)+ l) = f(i,:)*d(j,:,l)'; 
         elseif (m==i), 
           J(k,3*(m-1)+ l) = f(j,:)*d(i,:,l)';          
         end
       end
     end    
     k = k + 1;      
   end
end

% Note that the analytical Jacobian is close to the finite difference
% approximation 

end

