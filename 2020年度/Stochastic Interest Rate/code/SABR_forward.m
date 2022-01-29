%--------------------------------------------------------------------------
% SABR_forward - Simulate smile dymanics 
% SABR Model predicts smile/skew moves in the same direction as the value of the
% underlying; as the forward increases the Black implied volatility curve 
% moves to the right. 
%--------------------------------------------------------------------------

function SABR_forward

%f = 0.0652;
%beta = 0.5;
%sigma_0 = 0.02952;
%rho = -0.25;
%nu = 0.61;
%T = 0.50137;
f = 0.05;
beta = 0.5;
sigma_0 = 0.03;
rho = -0.25;
nu = 0.5;
T = 1;

% Need to add eps to forward TO avoid "Divide by zero" in SABR formula.
eps = 10^-8; 

%f = [0.0450    0.0470    0.0490    0.0510    0.0530    0.0550    0.050    0.0450] + eps;
f = [0.0450    0.0470    0.0490    0.0510    0.0530    0.0550] + eps;
f = [0.045 0.05] + eps;
Start = 0.04; 
End = 0.065;
  
for j = 1:length(f),
    
  K = [Start:0.001:End];
  for i = 1:length(K),
    sigma(i) = formula(f(j),K(i),sigma_0,beta,nu,rho,T*(1-j/252));        
  end
  set(gca,'Ydir','reverse')
  plot3(K,T*(1-j/252)*ones(1,length(K)),sigma,'k-');
  hold on; 

  % At-the-money implied volatilty
  K(i) = f(j);   
  ATM = formula(f(j)+eps,K(i),sigma_0,beta,nu,rho,T*(1-j/252));  
    
  % Path of ATM as the forward moves 
  backbone(j,:) = [f(j),T*(1-j/252),ATM];
end

%plot3(backbone(:,1),backbone(:,2),backbone(:,3),'--ko','LineWidth',1,...
%                'MarkerEdgeColor','k',...
%                'MarkerFaceColor','w',...
%                'MarkerSize',5)

plot3(backbone(:,1),backbone(:,2),backbone(:,3),'ko','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5)            
grid on;
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12')
xlabel(['Strike'],'FontSize',16,'Color','k');
ylabel(['Time To Expiry'],'FontSize',16,'Color','k'); 
zlabel(['Implied volatility'],'FontSize',16,'Color','k'); 
%title(['Smile Dynamics'],'FontSize',14,'Color','k');
view(45,25);
view(0,0); grid off; box on; 
end

function [sigma] = formula(f,K,sigma_0,beta,nu,rho,T) 
    A = sigma_0/(((f*K)^((1-beta)/2))*(1 + ((1-beta)^2)*(((log(f/K))^2)/24)) +  ((1-beta)^4)*(((log(f/K))^4)/1920));
    B = 1 + ((((1-beta)^2)/24)*(sigma_0^2)/((f*K)^(1-beta)) + (rho*beta*nu*sigma_0)/(4*(f*K)^((1-beta)/2)) + (2-3*rho^2)*(nu^2)/24)*T; 
    z = (nu/sigma_0)*((f*K)^((1-beta)/2))*log(f/K);
    chi = log((sqrt(1 - 2*rho*z + z^2) + z - rho)/(1-rho));
    sigma = A*(z/chi)*B;    
end 
