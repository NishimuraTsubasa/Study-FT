%--------------------------------------------------------------------------
% Daragh McInerney, 2013
% Use historical USD interest rate data to derive USD forward rates along with 
% a sample correlation matrix 
%--------------------------------------------------------------------------

clear all;
%--------------------------------------------------------------------------
% Use H15 data downloaded from the federal reserve web-site 
% H15 data: International Swaps and Derivatives Association (ISDA®) mid-market par swap rates. 
% Rates are for a Fixed Rate Payer in return for receiving three month LIBOR, and are based 
% on rates collected at 11:00 a.m. Eastern time by Garban Intercapital plc and published 
% on Reuters Page ISDAFIX
%--------------------------------------------------------------------------
filename = 'USD_IR_FRB_H15_Mar2013.csv';
fid = fopen(filename,'r');

% Load the data setting headerLines to 1 to skip the first line in the file 
C = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]','delimiter',',','headerLines',1);
dates = C{1}; 
dates = datenum(dates,'dd/mm/yyyy'); % convert dates to date-numbers using correct date format: 03/07/2000
% Note date-number representation used in Excel differs from that used in used in Matlab

data  = cell2mat(C(2:end)); 
fprintf(' Data in file %s is from %s to %s\n',filename, datestr(dates(1,:)),datestr(dates(end,:)));

% fseek(fileID, offset, origin) sets the file position indicator offset bytes from origin in the specified file.
fseek(fid, 0, 'bof'); 

% Data contains some columns with only NaN's. Find index of first column
% consisting only of NaN's. 
offset = length(find(isfinite(data(end,:)))) + 1;
% C = textscan(fid, 'format', N) reads data from the file, using the format N times, where N is a positive integer. 
C = textscan(fid, '%s', offset, 'delimiter', ',');

% transform expiries to integer number of years 
expiries = transpose(C{1}(2:offset));   
for i = 1:length(expiries),
  expiries{i} = str2num(expiries{i}(1:end-1));
end
expiries = cell2mat(expiries); 

if 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% liquidity effects may lead to a sharp drop in correlation between rates:
% 
% %usd_ir_data data contains swap rates for maturities [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 20 25 30 40 50]
% H15 data contains (mid-market par) swap rates for maturities [1 2 3 4 5 7 10 30] to ensure that 
% we only include "liquid" swap rates in the construction of the discount curve we select a subset of 
% the avialable rates: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set of liquid swap rates 
set = [1 2 3 4 5 7 10 12 15 20 30]; 
j = 1; 
for i = 1:length(set),
   index = find(expiries == set(i)); 
   if ~isempty(index), 
       indices(j) = index; 
       j = j + 1;
   end
end
data = data(:,indices); 
expiries = expiries(indices); 
end 

% Filter out non-trading days. Important for H15 data where non-trading weekdays/holidays are included.  
indices = find(sum(isnan(data),2) < length(expiries)); 

IRData.dates = dates(indices);
IRData.swapRates = 0.01*data(indices,:); 
IRData.expiries = expiries;

if 1 % Data check 
%--------------------------------------------------------------------------
% For data set: usd_ir_data
% Only after Aug 1998 do we have a series consisting of only trading days. 
% Suspect data? Even though non-trading days are included in the early part 
% of the time series (before Aug 1998) the data is not forward or backfilled 
% for those dates. 
%--------------------------------------------------------------------------

% Find the first instance in the series where non-trading days (that were not weekdays) were removed.
index = min(find(diff(IRData.dates) > 3));
fprintf(' Market closed from end of %s (%s) until beginning of %s (%s).\n', datestr(IRData.dates(index)), ... 
         datestr(IRData.dates(index),'ddd'),datestr(IRData.dates(index+1)),datestr(IRData.dates(index+1),'ddd')); 

% Date(s) where the data is identical to that on the previous date -
% possible indication of forward/backfilled data? 
indices = find(sum(abs(diff(IRData.swapRates(:,1:length(IRData.expiries)),1,1)),2) < 0.00001) + 1; %
for i = 1:length(indices), 
  fprintf(' Warning: Repeated/Identical data found on %s (%s) \n',datestr(IRData.dates(indices(i))), datestr(IRData.dates(indices(i)),'ddd'));
end

end

clear C dates data expiries 

%--------------------------------------------------------------------------
% Select start-date and end-date for forward rate time-series
%--------------------------------------------------------------------------

%date1 = '05-Aug-1987'; % start of usd_ir_data time-series
%date1 = '03-Jul-2000'; % start of H15 time-series 

startDate = '31-Jul-2004';
%startDate = '31-Jul-2008';

%endDate = '31-Jul-2009'; 
%endDate = '25-Mar-2011'; 
endDate = '15-Mar-2013';

startIndex = min(find(IRData.dates >= datenum(startDate)));
endIndex = min(find(IRData.dates >= datenum(endDate)));
fprintf(' Extracting rates from %s to %s \n',datestr(IRData.dates(startIndex)), datestr(IRData.dates(endIndex)));

%--------------------------------------------------------------------------
% Assume a set of dates: 0 < T_{0} < ... < T_{N} with an accural period of 1 year. 
% By default N = 20 this allows us to derive 21 forward rates, the first rate 
% expiring at T_{0} = 1 year and the 20th at T_{N-1}. In the notation used in 
% the text-book the rate expiring at time $T_{0}$; F(T_{0};T_{0},T_{1}) is denoted 
% by F_{1}(T_{0}), in other words the forward rate index refers to the index of 
% the date at the end of the accrual period. 
%--------------------------------------------------------------------------
accrual = 1; 
fprintf(' We fix the accrual period to be 1 year by default\n'); 
N = 21; 
fprintf(' Calculate discount factors out to %d years\n',N); 

swapRates = IRData.swapRates(startIndex:endIndex,:); 
dates = IRData.dates(startIndex:endIndex);
expiries = IRData.expiries; 

%-------------------------------------------------------------------------- 
% Transform to weekly data. Could use only weekly (Friday) data  
%--------------------------------------------------------------------------
if 1 
  % Take only weekly observations. Filter for all Friday data: 
  % indices    = strmatch('Fri',datestr(dates,'ddd'));
  % Or alternatively filter for market days just before week-ends
  indices    = find(diff(dates) >= 3);  
  % disp(datestr(dates(indices),'ddd'));
  dates      = dates(indices);
  swapRates = swapRates(indices,:); 
end

T = size(swapRates,1); 
fprintf(' Filter for trading days just before week-ends\n'); 
fprintf(' Use %d days in total from our sample data from %s to %s\n',T,datestr(dates(1)),datestr(dates(end))); 

%--------------------------------------------------------------------------
% Necessary to perform a cubic spline interpolation between spot rates 
% to estimate intermediate rates not present in the data. 
% Need to do this for H.15 data. 
%
% Rebonato & McKay p.149: In order to reduce oscillations away from the
% nodes induced by the cubic spline procedure it is also common to use
% natural splines (splines for which the second derivatives at both ends of
% the interpolation specture are set to zero)
%--------------------------------------------------------------------------

interpRates = NaN*zeros(length(swapRates),N);
for i = 1:length(swapRates),
  % Only use indices of swap rates that are valid (i.e. check for NaN or Inf)
  indices = find(isfinite(swapRates(i,:)));
   
  % interRates(i,:) = spline(IRData.expiries(indices),swapRates(i,indices),1:N);
  cs = spline1d(expiries(indices),swapRates(i,indices),1e30,1e30);
  interpRates(i,:) = ppval(cs,1:N); 
end

% Test: Plot interpolated rate-curve for a specific date
if 0
  date = input('Input date (08-Mar-2013): ','s');
  index = min(find(dates >= datenum(date)))
  
  indices = find(isfinite(swapRates(index,:)));
  plot(expiries(indices),swapRates(index,indices),'o',1:N,interpRates(index,:),'x-');
  fprintf(' Plot of interpolated rate-curve for %s\n',datestr(dates(index)));
  
  if 0
  %--------------------------------------------------------------------------  
  % Plot piecewise cubic spline function to test implementation 
  %--------------------------------------------------------------------------  
  %cs = spline(expiries(indices),swapRates(index,indices));
  cs = spline1d(expiries(indices),swapRates(index,indices),1e30,1e30);
  
  i = 6;    
  x = linspace(expiries(indices(i)),expiries(indices(i+1)),11);   
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
  fprintf('% 1.10f\n', 6*cs.coefs(end,1)*(expiries(indices(end))-expiries(indices(end-1))) + 2*cs.coefs(end,2));

  % Matlab implements "Not-A-Knot" end point condition - the third derivatives 
  % at the end points match
  if (abs(diff(cs.coefs(1:2,1))) < 1e-16) && abs((diff(cs.coefs(end-1:end,1))) < 1e-16) 
    fprintf(' "Not-A-Knot" constriant is satisfied\n');
  end  
  hold off
  pause;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Libor rates from set of swap rates by boot-strapping procedure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
discountFactors = zeros(T,N);
level = zeros(T,N); % PV01

discountFactors(:,1) = 1./(1 + accrual*interpRates(:,1));
level(:,1) = accrual*discountFactors(:,1); 

% Calculate discount factors 
for i = 2:N, 
  discountFactors(:,i) = (1 - interpRates(:,i) .* level(:,i-1)) ./ (accrual*interpRates(:,i) + 1); 
  level(:,i) = level(:,i-1) + accrual*discountFactors(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use log-differences or SABR-compatible changes in the forward rate?   
% 
% Note: For the time-series of 12-month Libor forward rates used to create the correlation matrix 
% the time-to-maturity is fixed (constant residual time to maturity) and the maturity 
% changes so at the end of the time series (5-Nov-2007) the first forward rate would 
% would still have 1-year to maturity.    
%
% Another approach (Brigo and Mercurio, chapter 6) is one where the maturity of each forward 
% rate is fixed and the time-to-maturity decreases while moving through time. Therefore, the 
% first forward rate in the family would expire in 1-year so only 1-year of data would be used. 
% This approach would seem to be preferable as forward rates in the LMM are characterized by 
% a fixed maturity as time passes. 
%
% In Rebonato & McKay (2009) and Longstaff, Santa-Clara and Schwartz (2001)  
% the time-series of Libor forward rates have a fixed time-to-maturity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Create a total fo N-1 forward rates for a set of N discount factors  
forwardRates = zeros(T,N-1);
for i = 2:N, 
  forwardRates(:,i-1) = (1/accrual)*((discountFactors(:,i-1) ./ discountFactors(:,i)) - 1); 
end
fprintf(' From our set of %d discount factors we create %d forward rates\n',N,N-1); 

% Plot Time-series of '1-year' to '10-year' forward rates
if 1
  % figure('Position',[10 80 800 400]) 
  figure('Position',[10 80 800 400]) 
  %set(0,'DefaultAxesFontName', 'Times New Roman')
  %set(0,'DefaultAxesFontSize', 12)
  
  colormap gray % 0.03 - 0.08 
  
  % for clarity skip some data points 
  skip = round(length(dates)/100);
  indices = 1:skip:length(dates);  
  
  for i = 1:10, 
      plot3(i*ones(1,length(indices)),dates(indices),forwardRates(indices,i),'Color',[0.2 0.2 0.2]);  
      hold on;
  end 
  xlabel(['Expiry'],'FontSize',16,'Color','k');
  zlabel(['Forward rate'],'FontSize',16,'Color','k'); 
  %title(['Time-Series of Forward Rates'],'FontSize',14,'Color','k');
  ylabel(['Date'],'FontSize',16,'Color','k');  
  %datetick('y',10); % 10 'yyyy'
  datetick('y',28,'keepticks'); % 28 'mmmyyyy'
  grid on; 
  box on;
  axis([1 10 dates(1) dates(end) 0.00 0.065])
  view(-90,0); 
    
end 

%--------------------------------------------------------------------------
% Choice of Beta has relatively limited impact on the correlation of 
% changes in Libor Forward rates  
%--------------------------------------------------------------------------
Beta = 1.0; % for Log-Normal model 
%Beta = 0.0; % Normal model  
%Beta = 0.5; % SABR (CEV = 0.5) model 

forwardRateDifferences = (forwardRates(2:end,:) - forwardRates(1:end-1,:)) ./ (forwardRates(1:end-1,:).^Beta);

% rho = corrcoef(X) returns a matrix rho of correlation coefficients calculated 
% from an input matrix X whose rows are observations and whose columns are variables
rho = corrcoef(forwardRateDifferences);

if 1
    figure
    colormap gray
    mesh(rho);
    axis([1 20 1 20 0.3 1])
    zlabel(['Correlation'],'FontSize',14,'Color','k'); 
    xlabel(['Forward-rate expiry'],'FontSize',14,'Color','k');
    ylabel(['Forward-rate expiry'],'FontSize',14,'Color','k');
    %title(['The Sample Correlation Matrix'],'FontSize',14,'Color','k');

    % Make up one column of the colormap. Now the colormap goes from 0 at index 1
    % corresponding to gray level 24, to 1 at index 256 corresponding to gray level 256.
    cmap1 = [linspace(0, 0.3, 256-0)];
    % Make up a gray colormap by copying this to all 3 columns.
    cmap3 = [cmap1; cmap1; cmap1]';
    % Apply the colormap to the image.
    colormap(cmap3); 
    view(-40,30)
end

% save variable rho in file rho_MMMYYYY.mat 
save rho_Mar2013.mat rho; 

return;

if 0
%--------------------------------------------------------------------------
% The correlation between the first and other forward rates between the second and 
% the others, between the third and the others and between the fourth and the others 
%--------------------------------------------------------------------------
plot(rho(:,1),'b-o'); % Circle
hold on;
plot(rho(:,2),'b-*'); % Asterix
plot(rho(:,3),'b-s'); % square
plot(rho(:,4),'b-d'); % diamond
 
ylabel(['Correlation '],'FontSize',12,'Color','k'); 
title(['Correlation of changes in Libor Forward Rates  (Figure 1)'],'FontSize',14,'Color','k');
end

if 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The correlation between the first five and subsequent forward rates: The first 
% series shows the correlation between the first forward rate (expiring in
% 1 year) and the forward rates expring in two, three,..., years time.
% The second series shows the correlation between the second forward rate (expiring in
% 2 years) and the forward rates expring in three, four,...,years time.

% Empirical features: 
% (i) The plot of the correlation between the front Libor forward and the later-expring ones is convex  
% (ii) for the later-expiring rates the correlation displays negative convexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off;
plot(rho(1:end-4,1),'b-o'); % Circle
hold on;
plot(rho(2:end-3,2),'b-*'); % Asterix
plot(rho(3:end-2,3),'b-s'); % square
%plot(rho(4:end-1,4),'b-d'); % diamond
%plot(rho(5:end,5),'b-^');% Upward-pointing triangle

ylabel(['Correlation '],'FontSize',12,'Color','k'); 
title(['Correlation of changes in Libor Forward Rates (Figure 2)'],'FontSize',14,'Color','k');
end

if 0    
[V,D] = eig(rho);

% Order the eigenvalues from highest to lowest 
[E,I] = sort(diag(D),1,'descend');

% Create associated matrix of eigenvectors (orthonormal set: U*U' = I)
U = V(:,I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Factors driving the term structure: 
% 
% The first factor generates (approximately) parallel shifts in the term structure. 
% The second factor generates shifts in the slope, while the third generates movements 
% in the term structure where both short-term and long-term rates move in an opposite 
% direction to mid-term rates. The fourth factor effects the shape of the short-end 
% of the term structure (Longstaff, Santa-Clara and Schwartz (2001)). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Percentage of variance explained by each eigenvalue 
for i = 1:5,
  fprintf('Eignevalue: %d Contribution %2.1f%%\n', i, 100*E(i)/length(E));
end
fprintf('First four factors account for %2.1f%% of the total variability in term structure.\n', 100*sum(E(1:4))/length(E)) ;

% Plot (normalised) eigenvectors associated with 4 largest eigenvalues 
figure;
hold on;
box on;
plot(U(:,1),'k-'); % parallel shift factor 
plot(U(:,2),'b-'); % slope 
plot(U(:,3),'r-'); % curvature
plot(U(:,4),'-', 'Color',[0.5 0.5 0.5]); % fourth factor 
xlabel(['Horizon'],'FontSize',12,'Color','k'); 
title(['The First Four Eigenvectors Of The Historical Correlation Matrix'],'FontSize',12,'Color','k');
end

% Notes: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eurodollars are deposits denominated in US dollars at banks outside the United States, 
% and thus are not under the jurisdiction of the Federal Reserve. Consequently, such deposits 
% are subject to much less regulation than similar deposits within the United States, allowing 
% for higher margins.
%
% For more than 30 years, the London interbank market has been a major source of liquidity for 
% banks seeking to fund U.S. dollar-denominated loans. The market has worked efficiently to 
% enable banks in need of liquidity to obtain deposits of U.S. dollars, either on an 
% overnight basis or for fixed terms (typically, one, two, three or six months), from other 
% banks with excess liquidity.
%
% Those deposits are sometimes referred to as “Eurodollar deposits”; and the interest rate 
% on Eurodollar deposits is referred to as “LIBOR” or the “LIBO Rate”, standing for the London 
% interbank offered rate.1 The British Bankers’ Association (BBA) oversees the public 
% quotation of LIBOR by designating a panel of 16 U.S. and non-U.S. banks2 to furnish their 
% rates to Reuters, which publishes the average after eliminating the highest and lowest rates.
%
% For the data in H.15 data set: 
% Instrument: "Eurodollar deposits (London)"
% Maturity: "6-month"
% Frequency: "Business day"
% Description: "U.S. -- SHORT-TERM INTEREST RATES: DAILY 6-MONTH EURO-DOLLAR DEPOSIT RATE"
% Note: "Annualized using a 360-day year or bank interest."
% Note: "Bid rates for Eurodollar deposits collected around 9:30 a.m. Eastern time."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CME Eurodollar futures prices are determined by the market’s forecast of 
% the 3-month USD LIBOR interest rate expected to prevail on the settlement date. 
% The settlement price of a contract is defined to be 100.00 minus the official 
% British Bankers Association fixing of 3-month LIBOR on the contract settlement date. 
% For example, if 3-month LIBOR sets at 5.00% on the contract settlement date, 
% the contract settles at a price of 95.00.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
