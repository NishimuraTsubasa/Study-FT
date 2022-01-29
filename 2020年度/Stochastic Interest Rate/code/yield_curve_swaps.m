%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daragh McInerney, 2013
%
% The construction of the discount/yield cure is not well posed in a mathematical sense. 
% We need to determine an infinite number of discount factors from a relatively small number 
% of market instruments. Therefore, we must use some form of interpolation. We make the 
% assumption that Instantaneous forward rates are constant between the benchmark maturities
% Here we just use swaps to construct the curve.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
% Interpolate swap-rates? 
INTERP_RATES = false; 
%startDate = '15-Mar-2013';
startDate = '18-May-2011';
%startDate = '25-Feb-2010';
%startDate = '18-Apr-2007';
%startDate = '13-Nov-2002'; 
%startDate = '15-Nov-2000'; 
%clear all;
format long;

if 1 
% load curve data from a csv file
%file_name = 'usd_curve_18Apr2007_ex1.csv';
file_name = 'usd_curve_18May2011_ex1.csv';
%file_name = 'usd_curve_25Feb2010.csv';
%file_name = 'usd_curve_15Jan2010.csv';
%file_name = 'usd_curve_13Nov2002.csv';
fid = fopen(file_name,'r');

% Set headerLines to 1 to skip the first line in the file 
C = textscan(fid,'%s%f%s%*[^\n]','delimiter',',','headerLines',1);

% Swap-rates 
j = 1;
swapRates = {}; 
for i = 1:length(C{1}) 
  if strcmpi(C{1}{i}(end),'y'), 
    if strcmpi(C{3}{i},'TRUE'),  
      swapRates(j,:) = {str2double(C{1}{i}(1:end-1)), 0.01*C{2}(i)}; 
      j = j + 1; 
    end 
  end   
end 
swapRates = cell2mat(swapRates); 
startDate = [file_name(11:12),'-',file_name(13:15),'-',file_name(16:19)];
fprintf('Date: %s (%s)\n',datestr(startDate), datestr(startDate,'ddd'));
end

if 0
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
dateNumbers = datenum(dates,'dd/mm/yyyy'); % convert dates to date-numbers using correct date format: 03/07/2000
% Note date-number representation used in Excel differs from that used in used in Matlab
data  = cell2mat(C(2:end)); 
fprintf(' Data in file %s is from %s to %s\n',filename, datestr(dateNumbers(1,:)),datestr(dateNumbers(end,:)));

% fseek(fileID, offset, origin) sets the file position indicator offset bytes from origin in the specified file.
fseek(fid, 0, 'bof'); 
% Data contains some columns with only NaN's. Find index of first column consisting only of NaN's. 
offset = length(find(isfinite(data(end,:)))) + 1;
% C = textscan(fid, 'format', N) reads data from the file, using the format N times, where N is a positive integer. 
C = textscan(fid, '%s', offset, 'delimiter', ',');

% transform expiries to integer number of years 
for i = 2:length(C{1}),
  swapRates(i-1,1) = str2double(C{1}{i}(1:end-1));
end

% extract swap rates for relevant start-date
index = find(dateNumbers == datenum(startDate));
swapRates(:,2) = data(index,1:offset-1) ./ 100; 
fprintf(' Extract %d swap rates for start-date: %s\n',length(swapRates),datestr(dateNumbers(index))); 

end

forwardRates = []; 
accrualDateNumbers = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boot-strapping the swap-rates: 
% Need to compute forward curve from our set of forward rates to estimate the discount factors 
% to use in the boot-strapping procedure for the first liquid swap rate
%
% We begin by interpolating the forward rates - there is no standard for interpolation 
% and various schemes have been proposed: 
% (i) Constant instantaneous forward rate: Instantaneous forward rates are assumed
% constant between the benchmark maturities ... 
% (ii) Linear instantaneous forward rate: Instantaneous forward rates are assumed
% linear between the benchmark maturities and continuous throughout. 
% (iii) Quadratic instantaneous forward rate: Instantaneous forward rates are assumed 
% quadratic between the benchmark maturities and continuously differentiable throughout. 
% This requires matching the values *and* the first derivatives of the instantaneous 
% forward rate at the benchmark maturities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spot starting swaps based on LIBOR start date 2 business days
% from the current date and mature and pay interest on anniversary dates that use 
% the same modified following business day conventions as the LIBOR index. Interest
% is usually computed on an act/360 day basis on the floating side of the swap and
% on 30/360 day basis in the fixed side. Typically, fixed payment dates
% (“coupon dates”) are semiannual (every 6 months), and floating payment dates are
% quarterly (every 3 months) to correspond to a 3 month LIBOR.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = datestr(busdate(busdate(startDate,1),1));
swapAccrualDates{1} = d;

% Total number of set/payment dates on the fixed side 
N = 2*swapRates(end,1) + 1;
for i = 2:N,           
    d = datestr(datevec(d) + [0 6 0 0 0 0]); 
    if isbusday(d), 
        swapAccrualDates{i} = d; 
    else % date falls on a weekday or holiday 
        if strcmpi(datestr(busdate(d,1),'mm'), datestr(d,'mm')), 
            swapAccrualDates{i} = datestr(busdate(d,1)); 
        else
            swapAccrualDates{i} = datestr(busdate(d,-1)); 
        end
    end         
end

disp(' '); 
for i = 1:N, 
    fprintf(' Coupon dates (semi-annual) for fixed side of swap: [%d] : %s\n', i-1, swapAccrualDates{i}); 
end

if 1
%--------------------------------------------------------------------------
% May be necessary to perform a cubic spline interpolation between spot rates 
% to estimate intermediate rates not present in the data. 
%
% Rebonato & McKay p.149: In order to reduce oscillations away from the
% nodes induced by the cubic spline procedure it is also common to use
% natural splines (splines for which the second derivatives at both ends of
% the interpolation specture are set to zero)
%--------------------------------------------------------------------------    
       
cs = spline1d(swapRates(:,1)',swapRates(:,2),1e30,1e30);
interpRates = ppval(cs,swapRates(1,1):swapRates(end,1));

interpRates = [swapRates(1,1):swapRates(end,1); interpRates]'; 

% Plot swap rates (circle) along with interpolated rates (cross) 
if 1
    figure;
    for i = 1:length(interpRates),
         index = find(interpRates(i,1) == swapRates(:,1));
         if isempty(index), 
             % plot((datenum(swapAccrualDates(2*interpRates(i,1) + 1)) - datenum(startDate))/365,interpRates(i,2),'bx-')
             plot((datenum(swapAccrualDates(2*interpRates(i,1) + 1)) - datenum(startDate))/365,interpRates(i,2),'kx-')
             hold on; 
         end
    end
    
    for i = 1:length(swapRates), 
        % plot((datenum(swapAccrualDates(2*swapRates(i,1) + 1)) - datenum(startDate))/365,swapRates(i,2),'ro')
        plot((datenum(swapAccrualDates(2*swapRates(i,1) + 1)) - datenum(startDate))/365,swapRates(i,2),'ko')
    end
    xlabel(['Maturity '],'FontSize',14,'Color','k');
    ylabel(['Swap rate'],'FontSize',14,'Color','k'); 
    %title(['Market Swap Rates (circle) Interpolated Rates (cross)'],'FontSize',14,'Color','k');
    xlim([0 interpRates(end,1) + 1])
    ylim([0 0.06])
end 
end

%--------------------------------------------------------------------------
% Extract forward rates from the set of Swap rates by a boot-strapping procedure
%--------------------------------------------------------------------------
% Array of discount Factors at each payment date on fixed side 
discountFactors = []; 

if INTERP_RATES,
    swapRates = interpRates;
end

tic; 
for i = 1:length(swapRates), 
   
    % Use all the pre-calculated forward rate data in boot-strapping procedure: 
    if ~isempty(accrualDateNumbers),
        startIndex = max(find(datenum(swapAccrualDates) <= accrualDateNumbers(end)));
        % fprintf('\n S_j: %s T_l: %s\n',swapAccrualDates{startIndex}, datestr(accrualDateNumbers(end)));
    else 
        % if accrualDateNumbers is empty then there are no pre-calculated forward rates so set startIndex to 0
        startIndex = 0; 
    end

    % Calculate *known* discount factors by integrating the forward Curve
    for j = (length(discountFactors) + 1):startIndex, 
      discountFactors(j) = exp(-sum((1/365)*forwardCurve(1:(datenum(swapAccrualDates(j)) - datenum(startDate))))); 
    end

    % denote swap rate by c 
    c = swapRates(i,2);
    % interest accures on a 30/360 day basis on the fixed side 
    tau = 0.5; 

    % Break the expresssion for the swap rate into a known (RHS) and unknown (LHS) part. 
    if startIndex > 0,
        M = (discountFactors(1) - c*tau*sum(discountFactors(2:end)))/exp(-sum((1/365)*forwardCurve(1:(accrualDateNumbers(end)-datenum(startDate))))); 
    else
        M = 1.0; 
        % set the first accrualDateNumber to the first reset date and set the startIndex to 1 
        accrualDateNumbers(end+1) = datenum(swapAccrualDates(1)); 
        startIndex = 1; 
    end

    % index of swapAccrualDates corresponding to swap rate maturity
    swapAccrualIndex = 2*swapRates(i,1) + 1;
    % time (mesaured in days from today) when swap rate matures 
    accrualDateNumbers(end+1) = datenum(swapAccrualDates(swapAccrualIndex)); 

    % Create unknown LHS and use Bisection method to compute the forward rate
    F = @(x) (1 + c*tau)*exp(-(accrualDateNumbers(end) - accrualDateNumbers(end-1))*x/365);

    for j = startIndex:swapAccrualIndex-2, 
        F = @(x) F(x) + c*tau*(exp(-(datenum(swapAccrualDates(j+1)) - accrualDateNumbers(end-1))*x/365)); 
        % fprintf(' %s to %s: index: %d\n',datestr(accrualDateNumbers(end-1)),swapAccrualDates{j+1},datenum(swapAccrualDates(j+1)) - accrualDateNumbers(end-1));     
    end

    forwardRates(end+1) = Bisection(M, 0.0, 1.0, 0.00001, F);

    % interpolate of the forward curve assuming *constant* forward rates  
    if length(forwardRates) > 1, 
        forwardCurve((accrualDateNumbers(end-1)+1 - datenum(startDate)):(accrualDateNumbers(end) - datenum(startDate))) = forwardRates(end);
    else
        forwardCurve(1:(accrualDateNumbers(end) - datenum(startDate))) = forwardRates(end);
    end

    % Linear interpolation of the forward rates 
    % xi = 1:t(end);
    % forwardCurve = interp1(t,forwardRates,xi,'linear');

    % piece-wise cubic spline interpolation 
    % xi = 1:t(end);
    % cs = spline1d(t,forwardRates,1e30,1e30);
    % forwardCurve = ppval(cs,xi);

end

clear discountFactors; 
if 1
disp(' ');     
for i = 1:length(swapRates), 
    fprintf(' %3d & %4.5f & %s & %2.7f \n', swapRates(i,1), swapRates(i,2), swapAccrualDates{2*swapRates(i,1)+1}, forwardRates(length(forwardRates) - length(swapRates) + i)); 
end 
disp(' '); 
end
toc;

if 1
    figure('Position',[10 300 600 400])

    % Compute discount curve from forward curve
    discountCurve = zeros(length(forwardCurve),1); 
    for i = 1:length(forwardCurve)
        discountCurve(i) = exp(-sum((1/365)*forwardCurve(1:i))); 
    end
    
    % compute the yield Curve or YTM from the discount Curve  
    yieldCurve = -log(discountCurve)' ./ ((1:length(discountCurve))/365); 
   
    hold on;
    %bar((1:t(end))/365,yieldCurve)
    % Change the color of the graph so that the bins are a shade of grey and the edges 
    % of the bins are also a shade of grey 
    %h = findobj(gca,'Type','patch');
    %set(h,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);   % 'Color',[0 0 1]
    %box on; 
    
    t = accrualDateNumbers - datenum(startDate);
    indices = 1:length(swapRates);
   
    % plot(t(indices)/365,forwardRates(indices),'ro')
    plot(t(indices)/365,forwardRates(indices),'ko','LineWidth',0.75)
    % plot((1:t(end))/365,forwardCurve,'b-')
    stairs((1:t(end))/365,forwardCurve,'k-')
    
    plot((1:t(end))/365,yieldCurve,'k-','LineWidth',0.75)
    
    %ylabel('Forward Rate Curve','FontSize',14,'Color','k'); 
    ylabel('Rate','FontSize',14,'Color','k'); 
    box on;
    hold on; 
    xlabel('Benchmark maturities','FontSize',14,'Color','k');     
    %title('The Instantaneous Forward Rate f(t,T)','FontSize',16,'Color','k'); 
    axis([ 0 30 0 0.06])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Libor rates from set of swap rates by boot-strapping procedure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

swapRates = interpRates; 

% accrual = swapRates(1,1); 
accrual = 1; 
fprintf(' We fix the accrual period to be 1 year by default\n'); 

level = 0.0; % PV01
discountFactors(1) = 1/(1 + accrual*swapRates(1,2));
level(1) = accrual*discountFactors(1); 

% Calculate discount factors 
for i = 2: length(swapRates), 
  accrual = swapRates(i,1) - swapRates(i-1,1);  
  discountFactors(i) = (1 - swapRates(i,2)*level(i-1))/(accrual*swapRates(i,2) + 1); 
  level(i) = level(i-1) + accrual*discountFactors(i);   
end

figure;
for i = 1:length(swapRates),
     plot((datenum(swapAccrualDates(2*swapRates(i,1) + 1)) - datenum(startDate))/365,discountFactors(i),'ko')
     hold on; 
end
plot((1:t(end))/365,discountCurve,'k-');    

ylabel('Discount curve','FontSize',14,'Color','k'); 
box on;
hold on; 
xlabel('Maturity','FontSize',14,'Color','k');     
%title('The Discount Curve B(t,T)','FontSize',16,'Color','k'); 
axis([ 0 30 0.2 1.0])
return;