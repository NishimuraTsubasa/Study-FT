%--------------------------------------------------------------------------
% Read implied volatility data and plot graphs 
%--------------------------------------------------------------------------
clear all;

filename = 'USDSwaptionVolsLN_18May2011.csv';
% use matlab function tblread to read implied swaption volatilities 
[data,varnames,casenames] = tblread(filename,','); 

tenors = []; 
for i = 1:size(varnames,1),
    index = strfind(varnames(i,:),'y');
    tenors(i) = str2num(varnames(i,1:index-1)); 
end

expiries = []; 
for i = 1:size(casenames,1),
    index = strfind(casenames(i,:),'y');
    expiries(i) = str2num(casenames(i,1:index-1)); 
end

vols = data; 

% Output the ATM vol matrix with Latex-type formmating 
fprintf('% s & ', ' '); 
for i = 1:length(tenors)-1, 
  fprintf('{\\bf %s} & ',strcat(num2str(tenors(i)),'y')); 
end
fprintf('{\\bf %s } \\\\ \n',strcat(num2str(tenors(i+1)),'y'));    
    
for i = 1:size(vols,1), 
    fprintf('{\\bf %s} & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f \\\\ \n',...
    strcat(num2str(expiries(i)),'y'), vols(i,:)); 
end 

%--------------------------------------------------------------------------
% Plot the ATM implied swaption vols 
%--------------------------------------------------------------------------
if 0 

hold on;
grid on;
% Note that 'x' corresponds to the columns of vols and 'y' corresponds to the rows.
mesh(tenors,expiries,vols)    

% Make up one column of the colormap. Now the colormap goes from 0 at index 1
% corresponding to gray level 24, to 1 at index 256 corresponding to gray level 256.
cmap1 = [linspace(0, 0.2, 256-0)];
% Make up a gray colormap by copying this to all 3 columns.
cmap3 = [cmap1; cmap1; cmap1]';
% Apply the colormap to the image.
colormap(cmap3); 

view(120,30);
xlabel(['Tenor'],'FontSize',14,'Color','k');
zlabel(['Implied volatility'],'FontSize',14,'Color','k'); 
title(['USD ATM Implied Swaption Volatilities'],'FontSize',14,'Color','k');
ylabel(['Expiry'],'FontSize',14,'Color','k');  
end 

%--------------------------------------------------------------------------
% Perform a cubic spline interpolation between ATM vols 
%
% To reduce oscillations away from the nodes induced by the cubic spline 
% procedure it is also a common to use natural splines (splines for which 
% the second derivatives at both ends of the interpolation specture are set 
% to zero)
%--------------------------------------------------------------------------

% new set of tenors 
interpTenors = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 12; 15; 20; 25; 30]; 
interpVols = NaN*zeros(length(expiries),length(interpTenors));
for i = 1:length(expiries),
  cs = spline1d(tenors,vols(i,:),1e30,1e30);
  interpVols(i,:) = ppval(cs,interpTenors); 
end

% plot interpolated vol-curve for a specific expiry
if 0
  
  expiry = input('Input expiry: ');
  index = min(find(expiries >= expiry))
  
  hold off;
  plot(tenors,vols(index,:),'o',interpTenors,interpVols(index,:),'x-');
  fprintf('Plot of interpolated vol curve for expiry: %d\n',expiries(index));
end

%--------------------------------------------------------------------------
if 1
% For a fixed swap maturity, interpolate along the option maturity axis: 
interpExpiries = [2; 3; 4; 5; 6; 7; 8; 9; 10; 12; 15; 20; 25; 30; 35; 40]; 
interpVolsExp = NaN*zeros(length(interpExpiries),length(interpVols));
for i = 1:length(interpVols),
  cs = spline1d(expiries,interpVols(:,i),1e30,1e30);
  interpVolsExp(:,i) = ppval(cs,interpExpiries); 
end
end


%--------------------------------------------------------------------------
% Plot the interpolated ATM implied swaption vols 
%--------------------------------------------------------------------------
if 1

hold off; 
mesh(interpTenors,interpExpiries,interpVolsExp);
grid on;

% Make up one column of the colormap. Now the colormap goes from 0 at index 1
% corresponding to gray level 24, to 1 at index 256 corresponding to gray level 256.
 cmap1 = [linspace(0, 0.2, 256-0)];
% Make up a gray colormap by copying this to all 3 columns.
cmap3 = [cmap1; cmap1; cmap1]';
% Apply the colormap to the image.
colormap(cmap3); 

view(120,30);
xlabel(['Tenor'],'FontSize',14,'Color','k');
zlabel(['Implied volatility'],'FontSize',14,'Color','k'); 
%title(['USD ATM Implied Swaption Volatilities'],'FontSize',14,'Color','k');
ylabel(['Expiry'],'FontSize',14,'Color','k');  
end

%--------------------------------------------------------------------------
% Plot the ATM implied  Cap vols 
%--------------------------------------------------------------------------
if 1 
figure
indices = [1 2 3 4 5 6 7 8 9 10 12 15 20]; 
caps = [10.38	15.11	16.04	16.19	16.18	16.12	15.97	15.76	15.55	15.36	14.90	14.38	13.80];
hold off; 
plot(indices,caps,'ko-')
xlabel(['Cap Maturity'],'FontSize',16,'Color','k');
%title(['USD ATM Implied Cap Volatilities'],'FontSize',22,'Color','k');
ylabel(['Implied Volatilites'],'FontSize',16,'Color','k'); 

figure; 
indices = [1 2 3 4 5 6 7 8 9 10 12 15 20]; 
caps = [90.30	82.67	66.97	55.31	46.90	41.42	37.68	35.03	33.13	31.64	29.45	27.11	24.95];
plot(indices,caps,'ko-')
xlabel(['Cap maturity'],'FontSize',16,'Color','k');
%title(['USD ATM Implied Cap Volatilities'],'FontSize',22,'Color','k');
ylabel(['Implied volatilites'],'FontSize',16,'Color','k');
end

%interpExpiries = 0.25:0.25:4 
%interpVols = NaN*zeros(length(interpExpiries));

%cs = spline1d(indices,caps,1e30,1e30);
%interpCaps = ppval(cs,Expiries); 
