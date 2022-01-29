function [cs, y2] = spline1d(x, y, yp1, ypn)
% Given arrays x[1...n] and y[1...n] containing a tabulated function y_i =
% f(x_i), with x_1 < x_2 < ... x_n, and given values yp1 and ypn for the
% first derivative of the interpolating function at points 1 and n,
% respectively, this routine returns an array y2[1...n] that contains the
% second derivatives of the interpolating function at the tabulated points
% x_i. If yp1 and/or ypn are equal to 10^30 or larger, the routine is
% signaled to set the corresponding boundary condition for a natural
% spline, with zero second derivative on that boundary. 

n = length(x); 
u = zeros(1,n);
y2 = zeros(1,n);

yp1 = 1e30; 
ypn = 1e30;

if yp1 > 0.99e30,% The lower boundary condition is set either to be "natural" or else 
    y2(1) = 0.0; % to have a specified first derivative.
    u(1) = 0.0;
else
    y2(1)  = -0.5; 
    u(1) = (3.0/(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1);
end

%
%
%
for i = 2:1:n-1, 
    sig = (x(i) - x(i-1))/(x(i+1) - x(i-1)); 
    p = sig*y2(i-1) + 2.0; 
    y2(i) = (sig - 1.0)/p; 
    u(i) = (y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1))/(x(i) - x(i-1));
    u(i) = (6.0*u(i)/(x(i+1) - x(i-1)) - sig*u(i-1))/p; 
end

if ypn > 0.99e30,   % The upper boundary condition is set either to be "natural" or else 
    qn = 0.0;       % to have a specified first derivative. 
    un = 0.0; 
else
    qn  = 0.5; 
    un = (3.0/(x(n) - x(n-1)))*(ypn - (y(n) - y(n-1))/(x(n) - x(n-1)));
end
y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0); 

for k = n-1:-1:1,
    y2(k) =  y2(k)*y2(k+1) + u(k); 
end

for i = 1:n-1, 
   coefs(i,1) = (y2(i+1) - y2(i))/(6*(x(i+1) - x(i))); 
   coefs(i,2) = y2(i)/2; 
   coefs(i,3) = (y(i+1) - y(i))/(x(i+1) - x(i)) - (x(i+1) - x(i))*y2(i)/2 -(x(i+1) - x(i))*(y2(i+1) - y2(i))/6;
   coefs(i,4) = y(i); 
end

% return the piecewise polynomial form of the cubic spline interpolant for 
% later use with ppval...
cs.form = 'pp';
cs.breaks = x; 
cs.coefs = coefs; 
cs.pieces = length(x) - 1;
cs.order = 4;
cs.dim = 1; 
end
    