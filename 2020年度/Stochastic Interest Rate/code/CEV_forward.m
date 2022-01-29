%--------------------------------------------------------------------------
% CEV_forward - Simulate smile dymanics of CEV model (or rather displaced-diffusion 
% approximation to CEV model. The CEV Model predicts smile/skew moves in the *opposite* 
% direction as the value of the underlying; as the forward increases the Black implied 
% volatility curve moves to the left.  
%--------------------------------------------------------------------------
clear all;
close all;
format long; 

% Parameters for the forward rate process 
F        = 0.05;
stepSize = 0.005;
f        = [F - stepSize F F + stepSize]; 
a        = 0.05; % diffusion_coefficient
r        = 0.05; 
T        = 1; 
strikes  = 0.03:stepSize:0.07;  

% Note: We set volatility to sigma*F/(F+a) to ensure that the price of the ATM
% call (or equivalently the implied vol) is the same for each choice of 
% displaced diffusion coefficient. 
sigma    = 0.2;
sigma_D  = sigma*F/(F + a); 
% Remark: If we set vol of the process 'sigma' equal to 'sigma*f(i)/(f(i)+a)'
% then the smile moves with the underlying ... 

for i = 2:length(f),  % In text-book we plot F = 0.05 and F = 0.055   
    for j = 1:length(strikes),
        % Calculate the call price using Black's formula for each strike. 
        % The forward and strike have been adjusted to account for the fact 
        % that the forward is modelled as a displaced diffusion process 
        px(i,j) = callbl(f(i)+a, strikes(j)+a, r, sigma_D,T);        
        % Use Bisection method to compute the implied volality given the price 
        % as calculated via the displaced diffusion process 
        Fn = @(x) callbl(f(i), strikes(j) , r, x, T);
        impliedVol(i,j) = Bisection(px(i,j), 0.01, 0.99, 1.0e-15, Fn);  
    end
    
    if i == 1,
       plot3(strikes,i*ones(1,length(strikes)),impliedVol(i,:),'k:');
    end
    
    if i == 2,
       plot3(strikes,i*ones(1,length(strikes)),impliedVol(i,:),'k-');
    end
    
    if i == 3,
       plot3(strikes,i*ones(1,length(strikes)),impliedVol(i,:),'k--');
    end    
    
    hold on;  
end

set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12')

%plot(strikes, impliedVol(1,:),'k-');
%hold on;
%plot(strikes, impliedVol(2,:),'k--');
%plot(strikes, impliedVol(3,:),'k:');

%grid on;
ylabel('Implied volatility','FontSize',16,'Color','k');
xlabel('Strike','FontSize',16,'Color','k');     
%title('Smile dynamics for Displaced Diffusion Process','FontSize',16,'Color','k'); 
%axis([strikes(1) strikes(end) impliedVol(end,end)-0.02 impliedVol(end,1)+0.02]); 
view(39,24);
view(0,0)
axis([0.04 0.06 2 3 0.18 0.22]);
zlabel('Implied volatility','FontSize',16,'Color','k');
box on;

% The change in the implied vol with respect to a change in the forward rate
d_impliedVol_d_F = (impliedVol(3,:) - impliedVol(1,:))/(2*(f(3) - f(1)));
% The change in the implied vol with respect to a change in the strike
d_impliedVol_d_K =  (impliedVol(2,3:end) - impliedVol(2,1:end-2))/(2*stepSize);

if 0 
figure
plot(d_impliedVol_d_F(2:end-1),'k-');
hold on; 
plot(d_impliedVol_d_K(1:end),'k--');

% relative magnitude 
c = d_impliedVol_d_K(1:end) ./ d_impliedVol_d_F(2:end-1);
end
