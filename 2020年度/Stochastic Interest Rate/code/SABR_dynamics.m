%--------------------------------------------------------------------------
% SABR_dynamics
% Model predicts smile/skew moves in the same direction as the value of the
% underlying; as the forward increases the Black implied volatility curve 
% moves to the right. 
%--------------------------------------------------------------------------
function SABR_dynamics 

f = 0.05;
beta =0.5; 
sigma_0 = 0.03;
%rho = -0.25; 
%nu = 0.5; 
T = 1.0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependence on simga_0
% The main effect of an increase in the initial volatility simga_0 is an almost parallel shift 
% upwards in the Black implied volatility smile. We also see  a mild steepening of the smile. 
%
% Dependence on beta 
% The slope of the smile steepens as beta decreses from 1 to 0 ... 
% We observe an incresase in out-of-the-money implied volatility and a decrease
% in in-the-money Black implied volatility. There is also a slight increase in curvature 
% while the smile tends to shift unpwards. 
% We set beta = 0.5. This is the market standard value. 
%
% Note: when varying beta need to adjust sigma_0 to ensure that ATM implied volatilites match 
% sigma_0 = sigma_0*f^(0.5)/f^beta
%
% Dependence on rho
% A decrease in rho from 0 to -0.5 causes the smile to steepen. 
%
% Dependence of nu 
% The curvature of the smile increases as v increases 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EPS = 1e-12;% introduce some small offset (EPS) to ensure that at no point in the loop below 
            % does f = K(i) this ensures that sigma(i) is well-defined   
Start = 0.040 + EPS; 
End = 0.065 + EPS;
K = [Start:0.001:End];

% change rho 
%if 1
rho = [0.0 -0.25 -0.5]; % vector of 3 values to test 
nu = 0.5;
for k=1:length(rho), 
    for i = 1:length(K),
        sigma(i,k) = formula(f,K(i),sigma_0,beta,nu,rho(k),T);
    end
end

set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12')

% plot first 3 test values 
plot(K,sigma(:,1),'k-');
hold on;
plot(K,sigma(:,2),'k--');
plot(K,sigma(:,3),'k:','LineWidth',2);

axis([Start End min(sigma(:,3))-0.01 max(sigma(:,3))+ 0.005]);
hold on; 
%axis([Start End 0.11 0.17])

xlabel(['Strike'],'FontSize',16,'Color','k');
ylabel(['Implied volatility SABR'],'FontSize',16,'Color','k'); 
%title(['Smile Dynamics as \rho Changes'],'FontSize',20,'Color','k');
%end

% change nu 
%if 0
nu = [0.2 0.35 0.5]; 
rho = -0.25;
figure
for k=1:length(nu), 
    for i = 1:length(K),
        % Note if z = 0 (ATM) then chi = 0 and therefore z/chi is not defined 
        sigma(i,k) = formula(f,K(i),sigma_0,beta,nu(k),rho,T); 
    end
end

set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12')

% plot first 3 test values 
plot(K,sigma(:,1),'k-');
hold on;
plot(K,sigma(:,2),'k--');
plot(K,sigma(:,3),'k:','LineWidth',2);

axis([Start End min(sigma(:,3))-0.01 max(sigma(:,3))+ 0.005]);
hold on; 
%axis([Start End 0.11 0.17])

xlabel(['Strike'],'FontSize',16,'Color','k');
ylabel(['Implied volatility SABR'],'FontSize',16,'Color','k'); 
%title(['Smile Dynamics as \nu Changes'],'FontSize',20,'Color','k');

end
%end 

function [sigma] = formula(f,K,sigma_0,beta,nu,rho,T) 
    A = sigma_0/(((f*K)^((1-beta)/2))*(1 + ((1-beta)^2)*(((log(f/K))^2)/24)) +  ((1-beta)^4)*(((log(f/K))^4)/1920));
    B = 1 + ((((1-beta)^2)/24)*(sigma_0^2)/((f*K)^(1-beta)) + (rho*beta*nu*sigma_0)/(4*(f*K)^((1-beta)/2)) + (2-3*rho^2)*(nu^2)/24)*T; 
    z = (nu/sigma_0)*((f*K)^((1-beta)/2))*log(f/K);
    chi = log((sqrt(1 - 2*rho*z + z^2) + z - rho)/(1-rho));
    sigma = A*(z/chi)*B;    
end 