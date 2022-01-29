%--------------------------------------------------------------------------
% Plot of Volatility smile/skew for displaced diffusion model  
%--------------------------------------------------------------------------
clear all; 
format long; 

F     = 0.05;
sigma = 0.2; 
T     = 1; 
% assume no discounting 
r     = 0.05; 

%--------------------------------------------------------------------------
% Note: Setting volatility to sigma*F/(F+a) ensures that the price of the ATM
% call (or equivalently the implied vol) is the same for each choice of 
% displaced diffusion coefficient 'a'. 
% Specifically for a = 0 (log-normal process)  
% callbl(F, F, r, sigma,T) should be appoximately equal to 
% callbl(F+a, F+a, r, sigma*F/(F+a),T)
%
% Function used by finanical tool-box: 
% px = blkprice(F+a, F+a, r, T, sigma*F/(F+a));
%--------------------------------------------------------------------------

diffusion_coeffs = [0.0 0.01 0.05 1.0]; 
% diffusion_coeffs = 0.05; 
strikes = [0.04:0.002:0.06];  

for i = 1:length(diffusion_coeffs),     
    for j = 1:length(strikes),

        % Calculate the call price using Black's formula for each strike. 
        % The forward and strike have been adjusted to account for the fact 
        % that the forward is modelled as a displaced diffusion process 
        px(i,j) = callbl(F + diffusion_coeffs(i), strikes(j) +diffusion_coeffs(i), r, sigma*F/(F+diffusion_coeffs(i)),T);

        % Use Bisection method to compute the implied volatility given the price 
        % as calculated via the displaced diffusion process 
        Fn = @(x) callbl(F, strikes(j) , r, x, T);
        impliedVol(i,j) = Bisection(px(i,j), 0.1, 0.5, 0.000001, Fn);  
    end
end

%--------------------------------------------------------------------------
% Calculate the implied vol under the normal model, i.e., we assume the forward
% follows a Gaussian process. We set volatility equal to sigma*F so that 
% the ATM price as calculated by the normal model matches the displaced 
% diffusion ATM price. 
%--------------------------------------------------------------------------

% For a displaced diffusion coefficent >= 1 the displaced diffusion process can 
% be appoximated by a Gaussian process. For a value >> 1 the appoximation is 
% almost exact. 
i = i + 1;
for j = 1:length(strikes),
    % The call price assuming a normal model for the forward
    px(i,j) = callbln(F, strikes(j), r, sigma*F,T);
    Fn = @(x) callbl(F, strikes(j) , r, x, T);
    impliedVol(i,j) = Bisection(px(i,j), 0.001, 1.0, 0.0000001, Fn); 
end

% Plot of implied volatilities 
for i = 1:size(impliedVol,1)-1, 
   plot(strikes, impliedVol(i,:),'k-');    
   hold on;
end

set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12')

plot(strikes, impliedVol(end,:),'k--');
text(0.0405,0.1980,'$\alpha = 0.0$','Interpreter','latex','FontSize',14)
text(0.0405,0.2012,'$\alpha = 0.05$','Interpreter','latex','FontSize',14)
text(0.0405,0.2060,'$\alpha = 0.5$','Interpreter','latex','FontSize',14)
text(0.0405,0.2125,'$\alpha = 1.0$','Interpreter','latex','FontSize',14)

ylabel('Implied volatility','FontSize',16,'Color','k');
box on;
xlabel('Strike','FontSize',16,'Color','k');     
% title('Volatility Skew for Displaced Diffusion Process','FontSize',16,'Color','k'); 
axis([strikes(1) strikes(end) impliedVol(end,end) impliedVol(end,1)]); 
