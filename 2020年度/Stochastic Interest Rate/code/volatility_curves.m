%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot instantaneous volatility curves for forward rates for different
% parameterisations in both 'normal' and 'stressed' market conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function volatility_curves

delta_t = 0.05; 
expiry  = [0.0:delta_t:10];
   
%a       = 0.01344948;   
%b       = 0.19085664;
%c       = 0.97462314;
%d       = 0.08089168;  

%a = -0.02;
%b = 0.3;
%c = 2.0;
%d = 0.14;

% series 4
%a  = -0.05
%b = 0.7
%c = 1.5
%d = 0.1 

a  = -0.1;
b = 0.5;
c = 1;
d = 0.1;
fprintf(' The maximum of the hump: %f\n',1/c - a/b); 

for i = 1:length(expiry),      
  volatility(i)   = (a + b*expiry(i))*exp(-c*expiry(i)) + d;
end
plot(expiry,volatility,'k-')
hold on;

%--------------------------------------------------------------------------
%  Plot curve when in stressed or 'excited' market conditions 
%--------------------------------------------------------------------------

% series 5
% a = 0.3
% b = 1.5
% c = 5
% d = 0.15

a = 0.25;
b = 1;
c = 5;
d = 0.15;

fprintf(' The maximum of the hump: %f\n',1/c - a/b); 

for i = 1:length(expiry),      
  volatility(i)   = (a + b*expiry(i))*exp(-c*expiry(i)) + d;
end
plot(expiry,volatility,'k--')

%axis([0.0 10.0 0.10 0.20])
ylabel(['Volatility'],'FontSize',12,'Color','k'); 
xlabel(['Time to expiry'],'FontSize',12,'Color','k');
%title(['The Instantaneous Volatility Curve'],'FontSize',12,'Color','k');

end

