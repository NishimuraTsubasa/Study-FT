%--------------------------------------------------------------------------
% Read implied volatility data and plot graphs 
%--------------------------------------------------------------------------
clear all;

filenames = 'USDvolsLN.csv';
fid = fopen(file_name,'r');
% Set headerLines to 1 to skip the first line in the file 

% tenors 
C = textscan(fid,'%s%s%s%s%s%s%s%s%s%s*[^\n]','delimiter',',','headerLines',1);
C = transpose(C);   
% transform tenors to integer number of years 
for i = 1:length(C)-1,
  tenors{i} = str2num(C{i+1}{1}(1:end-1));
end
tenors = cell2mat(tenors); 

i = 1; 
while (1), 
    C = textscan(fid,'%s%f%f%f%f%f%f%f%f%f*[^\n]','delimiter',',','headerLines',1);    
    if isnan(C{3})
        break;
    end
    % transform expiries to integer number of years 
    expiries(i) = str2num(C{1}{1}(1:end-1));
    for j = 1:length(C)-1,
        vols(i,j) = C{j+1}; 
    end
    i = i + 1; 
end

% output the ATM vol matrix with Latex-type formating 
fprintf('% s & ', ' '); 
for i = 1:length(tenors)-1, 
  fprintf('{\\bf %s} & ',strcat(num2str(tenors(i)),'y')); 
end
fprintf('{\\bf %s } \\\\ \n',strcat(num2str(tenors(i+1)),'y'));    
    
for i = 1:size(vols,1), 
    fprintf('{\\bf %s} & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f \\\\ \n',...
    strcat(num2str(expiries(i)),'y'), vols(i,:)); 
end 

% plot the ATM implied swaption vols 

% Make up one column of the colormap. Now the colormap goes from 0 at index 1
% corresponding to gray level 24, to 1 at index 256 corresponding to gray level 256.
cmap1 = [linspace(0, 0.2, 256-0)];
% Make up a gray colormap by copying this to all 3 columns.
cmap3 = [cmap1; cmap1; cmap1]';
% Apply the colormap to the image.
colormap(cmap3); 
hold on;
grid on;

% Note that 'x' corresponds to the columns of vols and 'y' corresponds to the rows.
mesh(tenors,expiries,vols)

view(120,30);
xlabel(['Tenor'],'FontSize',14,'Color','k');
zlabel(['Implied Volatility'],'FontSize',14,'Color','k'); 
title(['USD ATM Implied Swaption Volatilities'],'FontSize',14,'Color','k');
ylabel(['Expiry'],'FontSize',14,'Color','k');  
  
%--------------------------------------------------------------------------
%  Perform a cubic spline interpolation between ATM vols 
%
% Rebonato & McKay p.149: In order to reduce oscillations away from the
% nodes induced by the cubic spline procedure it is also a common to use
% natural splines (splines for which the second derivatives at both ends of
% the interpolation specture are set to zero)
%--------------------------------------------------------------------------

indices = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 12; 15; 20; 25; 30]; 
interp_vols = NaN*zeros(length(expiries),length(indices));
for i = 1:length(expiries),
% interp_vols(i,:) = spline(tenors,vols(i,:),1:N);
  cs = spline1d(tenors,vols(i,:),1e30,1e30);
  interp_vols(i,:) = ppval(cs,indices); 
end

% plot interpolated vol-curve for a specific expiry
if 1
  
  expiry = input('Input expiry: ');
  index = min(find(expiries >= expiry))
  
  hold off;
  plot(tenors,vols(index,:),'o',indices,interp_vols(index,:),'-');
  fprintf('Plot of interpolated vol curve for expiry: %d\n',expiries(index));
  
  if 0
  %--------------------------------------------------------------------------  
  % Plot piecewise cubic spline function to test implementation 
  %--------------------------------------------------------------------------  
  %cs = spline(tenors,vols(index,:));
  cs = spline1d(tenors, vols(index,:),1e30,1e30);
  
  i = 4;    
  x = linspace(tenors(i),tenors(i+1),11);   
  hold on;
  plot(x,ppval(cs,x),'r-');
  
  a = cs.coefs(i,:);
  % R = @(t) a(1)*(t-x(1))^3 + a(2)*(t-x(1))^2 + a(3)*(t-x(1)) + a(4);
  plot(x,a(1)*(x-x(1)).^3 + a(2)*(x-x(1)).^2 + a(3)*(x-x(1)) + a(4),'r*'); 

  % Do first and second derivative constriants match adjacent splines? 
  fprintf('% 1.10f\n', 3*a(1)*(x(end)-x(1))^2 + 2*a(2)*(x(end)-x(1)) + (a(3) - cs.coefs(i+1,3)));
  fprintf('% 1.10f\n', 6*a(1)*(x(end)-x(1)) + 2*(a(2) - cs.coefs(i+1,2)));

  % Are there end point constraints that set the (*second*) derivative to zero at both ends?
  fprintf('% 1.10f\n', cs.coefs(1,2));
  fprintf('% 1.10f\n', 6*cs.coefs(end,1)*(tenors(end)- tenors(end-1)) + 2*cs.coefs(end,2));

  % Matlab implements "Not-A-Knot" end point condition - the third derivatives 
  % at the end points match
  if (abs(diff(cs.coefs(1:2,1))) < 1e-16) && abs((diff(cs.coefs(end-1:end,1))) < 1e-16) 
    fprintf(' "Not-A-Knot" constriant is satisfied\n');
  end  
  hold off
  pause;
  end
end

%MattF Wilmott: 
%For a fixed swap maturity, ie interpolating along the option maturity axis, 
%I think simply linearly interpolating on \sigma^2 T  is reasonable.

if 1
indices2 = [2; 3; 4 5; 6; 7; 8;	9; 10; 12; 15; 20; 25; 30; 35; 40]; 
interp_vols2 = NaN*zeros(length(indices2),length(interp_vols));
for i = 1:length(interp_vols),
% interp_vols2(:,i) = spline(expiries,interp_vols(:,i),1:N);
  cs = spline1d(expiries,interp_vols(:,i),1e30,1e30);
  interp_vols2(:,i) = ppval(cs,indices2); 
end
end

if 1
  
% plot the ATM implied swaption vols 
hold off; 
mesh(indices,indices2,interp_vols2);
hold on;
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
zlabel(['Implied Volatility'],'FontSize',14,'Color','k'); 
title(['USD ATM Implied Swaption Volatilities'],'FontSize',14,'Color','k');
ylabel(['Expiry'],'FontSize',14,'Color','k');  
end


if 1 
    
indices = [1 2 3 4 5 6 7 8 9 10 12 15 20]; 
caps = [10.38	15.11	16.04	16.19	16.18	16.12	15.97	15.76	15.55	15.36	14.90	14.38	13.80];
hold off; 
plot(indices,caps,'ko-')
xlabel(['Cap Maturity'],'FontSize',20,'Color','k');
title(['USD ATM Implied Cap Volatilities'],'FontSize',22,'Color','k');
ylabel(['Implied Volatilites'],'FontSize',20,'Color','k'); 

figure; 
indices = [1 2 3 4 5 6 7 8 9 10 12 15 20]; 
caps = [90.30	82.67	66.97	55.31	46.90	41.42	37.68	35.03	33.13	31.64	29.45	27.11	24.95];
plot(indices,caps,'ko-')
xlabel(['Cap Maturity'],'FontSize',20,'Color','k');
title(['USD ATM Implied Cap Volatilities'],'FontSize',22,'Color','k');
ylabel(['Implied Volatilites'],'FontSize',20,'Color','k');



end

