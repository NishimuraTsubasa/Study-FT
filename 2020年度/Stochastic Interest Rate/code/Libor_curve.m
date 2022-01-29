%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daragh McInerney, 2011. 
%
% The construction of the discount/yield cure is not well posed in a mathematical sense. 
% We need to determine an infinite number of discount factors from a relatively small number 
% of market instruments. Therefore, we must use some form of interpolation. We make the 
% assumption that Instantaneous forward rates are constant between the benchmark maturities
% We use Deposit rates, Eurodollar futures and Swaps to construct the curve.
% 
% Preference should always be given to more liquid instruments (i.e. those with 
% a tighter bid/ask spread). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;
format long;
IS_VASICEK = false % to run vasicek calibration

% load curve data from a csv file
%file_name = 'usd_curve_18Apr2007.csv';
file_name = 'usd_curve_18May2011_ex2.csv';
%file_name = 'usd_curve_25Feb2010.csv';
%file_name = 'usd_curve_15Jan2010.csv';
%file_name = 'usd_curve_13Nov2002.csv';
fid = fopen(file_name,'r');

% Set headerLines to 1 to skip the first line in the file 
C = textscan(fid,'%s%f%s%*[^\n]','delimiter',',','headerLines',1);

% LIBOR deposits 
j = 1;
deposits = {}; 
for i = 1:length(C{1}) 
  if length(C{1}{i}) < 8,
    if strcmpi(C{3}{i},'TRUE'),     
      deposits(j,:) = {C{1}{i}, C{2}(i)};      
      j = j + 1;
    end
  else 
    break; 
  end   
end

% Eurodollar Futures 
start = i; 
j = 1;
futures = {};
for i = start:length(C{1}) 
  if length(C{1}{i}) > 8, 
    if strcmpi(C{3}{i},'TRUE'),  
      futures(j,:) = {C{1}{i}, C{2}(i)}; 
      j = j + 1;
    end
  else 
    break; 
  end   
end

% Swap-rates 
start = i; 
j = 1;
swapRates = {}; 
for i = start:length(C{1}) 
  if length(C{1}{i}) < 8, 
    if strcmpi(C{3}{i},'TRUE'),  
      swapRates(j,:) = {str2double(C{1}{i}(1:end-1)), 0.01*C{2}(i)}; 
      j = j + 1; 
    end
  else 
    break; 
  end   
end 

startDate = [file_name(11:12),'-',file_name(13:15),'-',file_name(16:19)];
fprintf('Date: %s (%s)\n',datestr(startDate), datestr(startDate,'ddd'));
forwardRates = []; 
accrualDateNumbers = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
% Deposit Rates: 
% Deposits pay simple interest using an Act/360 day basis and a
% modified following convenction (if a date falls on a weekend or holiday,
% then it is moved to the following business day. If this places the date
% into the next month, then it is moved to the previous business day). 
% Deposits settle two business days after the trade date. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settlement_day = datestr(busdate(busdate(startDate,1),1)); 
despositAccrualDates = cell(length(deposits),2); 
for i = 1:length(deposits),           
  tenor = deposits{i,1};          

  if strcmpi(tenor,'O/N')  
      despositAccrualDates(i,:) = {startDate,datestr(busdate(startDate,1))};
  else % must be monthly deposit 
      d = datestr(datevec(settlement_day) + [0 str2double(tenor(1)) 0 0 0 0]); 
      if isbusday(d), 
         despositAccrualDates(i,:) = {settlement_day, d};
      else % date falls on a weekday or holiday 
          if strcmpi(datestr(busdate(d,1),'mm'), datestr(d,'mm')), 
              despositAccrualDates(i,:) = {settlement_day, datestr(busdate(d,1))}; 
          else
              despositAccrualDates(i,:) = {settlement_day, datestr(busdate(d,-1))}; 
          end
      end         
  end
end

%--------------------------------------------------------------------------
% Discount Factors 
% 'O/N' rate: P(t,t_1), '1M' rate: P(t',t_2), ... '6M' rate: P(t',...) 
%--------------------------------------------------------------------------                   
for i = 1:length(deposits),  
   accrualDateNumbers(i) = datenum(despositAccrualDates{i,2});
   discountFactors(i) = 1/(1 + 0.01*deposits{i,2}*((datenum(despositAccrualDates{i,2}) - datenum(despositAccrualDates{i,1}))/360));
end

if ~isempty(deposits), 
    % 'O/N': f(t_1)
    forwardRates(1) = -log(discountFactors(1))/((datenum(despositAccrualDates{1,2}) - datenum(despositAccrualDates{1,1}))/365);

    % First monthly-deposit: f(t_2)
    forwardRates(2) = -log(discountFactors(2))/((datenum(despositAccrualDates{2,2}) - datenum(despositAccrualDates{2,1}))/365);
end

for i = 3:length(deposits), 
  forwardRates(i) = -(1/((1/365)*(accrualDateNumbers(i) - accrualDateNumbers(i-1))))*log(discountFactors(i)/discountFactors(i-1));
end
clear discountFactors 

if 1 
for i = 1:length(deposits), 
    fprintf('%4s & %2.5f & %s to %s & %2.7f\n', deposits{i,1}, deposits{i,2}, despositAccrualDates{i,1}, despositAccrualDates{i,2}, forwardRates(i)); 
end    
end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eurodollar Futures 
% The quoted price of a Eurodollar future corresponds to 100(1 - Rate)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Bootstrapping EONIA or 3m Euribor? 
% Euribor: large credit risk premium embedded in the rate for unsecured lending.
% Theoretical discount factors are based on the price of a riskless bond. EURIBOR rates, being by 
% definition offer rates, carry the credit risk premium. This premium is not very high (roughtly 
% corresponding to AA-rated counterparty), but still the risk on a 3 month loan is lower than that 
% on a 6 month loan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------        
%[BeginDates, EndDates] = thirdwednesday(Month, Year) computes the beginning and 
% end period date for a LIBOR contract (third Wednesdays of delivery months).
% BeginDates is the beginning of three-month period contract as specified by Month and Year.
% EndDates is the end of three-month period contract as specified by Month and Year.
%--------------------------------------------------------------------------

if ~isempty(futures), 
    [Y, M, D] = datevec(datenum(futures{1,1}));
    futuresAccruals = zeros(length(futures),2); 
    for i = 1:length(futures), 
        [futuresAccruals(i,1), futuresAccruals(i,2)] = thirdwednesday(M,Y);
        [Y, M, D] = datevec(datestr(futuresAccruals(i,2)));
    end

    % transform IMMdata to 3-month futures rate
    futuresRate = zeros(length(futures),1);
    for i = 1:length(futures)
      futuresRate(i) = (1 - futures{i,2}/100);   
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eurodollar futures (convexity-adjusted?)
%
% Need to get the implied forward rate from an interest rate futures price: 
% The Eurodollar futures price not only incorporates information on the relevant forward 
% rate (the input necessary for building an interest rate term curve), but also on what is 
% known as a convexity adjustment. 

% FRA settle at maturity (i.e. there are no intermediate cash flows) while the Euro-dollar 
% futures contract is marked to market daily by the exchange. For a long position in a futures 
% contract if rates fall profits will be reinvested at a lower rate while if rates rise 
% losses will need to be financed at a higher rate.
% 
% Lesniewski (Convexity): Daily cash flows in and out of the margin account. The investor’s 
% P&L is negatively correlated with the dynamics of interest rates: If rates go up, the price 
% of the contract goes down, and the investor needs to add money into the margin account, 
% rather than investing it at higher rates (opportunity loss for the investor). If rates 
% go down, the contract’s price goes up, and the investor withdraws money out of the margin 
% account and invests at a lower rate (again an opportunity loss for the investor).
%
% To compenstate for this daily mark to market effect the price of the futures contract is 
% less then the equivalent LIBOR forward. Therefore, the futures rate calculated from the 
% price of a Eurodollar futures contract is higher than the corresponding LIBOR forward. 
%
% The convexity adjustment depends on how one models the future evolution of interest rates. 
% Let f(t,T,S) to be the forward rate, at time t, for the period [T,S]. The convexity adjustment 
% is defined as the difference: E[f(T,T,S)] - f(t,T,S) where the expectation is taken under 
% the risk neutral measure. 
% Note that the expectation of the spot rate at time T is simply the furtures rate at 
% time t for the period [T,S]. 
%
% Apply convexity-Adjustment to Futures rate using the Hull-White short-rate model
%    dr = (theta(t) - ar )dt + sigma dW
% a is the mena-reversion rate and sigma is the standard deviation of the
% change in the short-rate expressed annually 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0 
% Default parameters for the Hull-White short rate model used by liborfloat2fixed 
a = 0.05;    
sigma = 0.015; 
%--------------------------------------------------------------------------
% Uri Ron: A Practical Guide to Swap Curve Construction: Bloomberg assumes a constant 
% mean reversion rate of 0.03 %...
%--------------------------------------------------------------------------

for i = 1:length(futures),  
    %----------------------------------------------------------------------
    % In function liborfloat2fixed:  
    % t1 = yearfrac(StartDate,ReqminEndDate,2);
    % 
    % This is not correct, t1 should be set equal to: 
    % t1 = yearfrac(StartDate,ReqminBgnDate,2);
    %----------------------------------------------------------------------
    t1 = (futuresAccruals(i,1) - datenum(startDate))/360;
    B0_t1   = (1 - exp(-a*t1))/a;

    PeriodLength = (futuresAccruals(i,2) - futuresAccruals(i,1))/360;
    Bt1_t2  = (1 - exp(-a*PeriodLength))/a;
    
   % The convexity adjustment 
    FwdAdjust = Bt1_t2 ./ (PeriodLength) .* ...
        (Bt1_t2 .* (1 - exp(-2*a*t1)) + 2*a*B0_t1.^2)*(sigma^2)/(4*a);
    
   % The futures rate is quarterly compounded so we need to convert to a continuosuly 
   % compounded rate before making the adjustment 
   futuresRate(i) = 4*log(1 + futuresRate(i)/4) - FwdAdjust;
   
   % The convexity adjusted forward is now converted to a quarterly compounded rate
   futuresRate(i) = 4*(exp(futuresRate(i)/4) - 1);
  
   fprintf(' The convexity adjustment for IMM future beginning at: %s is %2.6f\n', datestr(futuresAccruals(i,1)),(1 - futures{i,2}/100)-futuresRate(i));
end
end

start = length(deposits); 
for i = 1:length(futures), 
    % fprintf('Determine the curve out to: %s\n',datestr(futuresAccruals(i,2)));
    accrualDateNumbers(start+i) = futuresAccruals(i,2);
    % for i > 1: accrualDateNumbers(5+i-1) - futuresAccruals(i,1)  is equivalent to: 
    % futuresAccruals(i-1,2) - futuresAccruals(i,1) and is typically small 
        d1 = futuresAccruals(i,2) - futuresAccruals(i,1);
        d2 = accrualDateNumbers(start+i-1) - futuresAccruals(i,1);
        d3 = accrualDateNumbers(start+i) - accrualDateNumbers(start+i-1);   
    if accrualDateNumbers(start+i-1) - futuresAccruals(i,1) > 0, 
        %fprintf(' d2 > 0:  d1 = %d, d2 = %d, d3 = %d\n',d1, d2, d3); 
        %fprintf(' Futures rate:%2.7f\n',futuresRate(i))
        %fprintf(' Forward rate:%2.7f\n',forwardRates(start+i-1))
        forwardRates(start+i) = (log(1 + futuresRate(i)*((futuresAccruals(i,2) - futuresAccruals(i,1))/360)) - forwardRates(start+i-1)*((accrualDateNumbers(start+i-1) - futuresAccruals(i,1))/365))/((accrualDateNumbers(start+i) - accrualDateNumbers(start+i-1))/365); 
    else 
        %fprintf(' d2 <= 0:  d1 = %d, d2 = %d, d3 = %d\n',d1, d2, d3); 
        %fprintf(' Futures rate:%2.7f\n',futuresRate(i))
        forwardRates(start+i) = log(1 + futuresRate(i)*((futuresAccruals(i,2) - futuresAccruals(i,1))/360))/((futuresAccruals(i,2) - futuresAccruals(i,1))/365);
    end
    %fprintf(' forward-rates: %2.7f\n',forwardRates(start+i));
end

if 1 
disp(' ');     
for i = 1:length(futures), 
    fprintf(' %s & %2.5f & %s to %s & %2.7f & %2.7f \n', datestr(datenum(futures{i,1}),'mmm yyyy'), futures{i,2}, datestr(futuresAccruals(i,1)), datestr(futuresAccruals(i,2)), forwardRates(start+i),(1 - futures{i,2}/100)-futuresRate(i)); 
end    
end

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

% Constant forward rates 
if ~isempty(accrualDateNumbers), 
    forwardCurve(1:(accrualDateNumbers(1)-datenum(startDate))) = forwardRates(1);
    for i = 2:length(accrualDateNumbers),
        forwardCurve((accrualDateNumbers(i-1)+1-datenum(startDate)):(accrualDateNumbers(i)-datenum(startDate))) = forwardRates(i); 
    end
end
       
% Linear interpolation of the forward rates 
% xi = 1:t(end);
% forwardCurve = interp1(t,forwardRates,xi,'linear');

% piece-wise cubic spline interpolation 
% xi = 1:t(end);
% cs = spline1d(t,forwardRates,1e30,1e30);
% forwardCurve = ppval(cs,xi); 

% Plot forward rates derived from Deposit rates and futures 
if 0
    t = accrualDateNumbers - datenum(startDate);
    plot(t/365,forwardRates,'rx')
    hold on; 
    plot((1:t(end))/365,forwardCurve,'b-')
    xlabel('Benchmark maturities (red)','FontSize',12,'Color','k');
    ylabel('Forward rate (instantaneous)','FontSize',12,'Color','k'); 
    title('The instantaneous forward rate f(0,t)','FontSize',14,'Color','k');
end 

% Swap rate data in Jeff Greco: "Using a Bootstrap Procedure to Build a LIBOR Curve" 
% dose not match other market data (4Y: 3.161, 5Y: 3.491, 7Y: 4.033, 10Y: 4.555 )
if 0            
% usd_ir_data for 13-Nov-2002          
swapRates = { '1Y', 1.533396363;'2Y', 2.107579164;'3Y', 2.616492619;'4Y', 3.00591036;'5Y', 3.327838178;'6Y', 3.596671594;'7Y', 3.829552436; '8Y', 4.027433277;'9Y', 4.192814119;'10Y', 4.33819496 }; 
% FRB - H15 for 13-Nov-2002 (mid-market par swap rates)
swapRates = { '1Y', 1.55;       '2Y', 2.1;        '3Y', 2.61;       '4Y', 3.0;       '5Y', 3.32;                         '7Y', 3.84;                                            '10Y', 4.36 }; 
end
swapRates = cell2mat(swapRates); 

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

% Plot swap rates (red) along with interpolated rates (blue) 
if 0 
    figure;
    for i = 1:length(interp_rates),
         plot((datenum(swapAccrualDates(2*interpRates(i,1) + 1)) - datenum(startDate))/365,interpRates(i,2),'bx-')
         hold on; 
    end
    
    for i = 1:length(swapRates), 
        plot((datenum(swapAccrualDates(2*swapRates(i,1) + 1)) - datenum(startDate))/365,swapRates(i,2),'ro')
    end
    xlabel(['Time to expiry '],'FontSize',12,'Color','k');
    ylabel(['Swap rates'],'FontSize',12,'Color','k'); 
    title(['Market Swap rates (red) interpolated rates (blue)'],'FontSize',14,'Color','k');
    xlim([0 interpRates(end,1) + 1])
end 
end

%--------------------------------------------------------------------------
% Extract forward rates from the set of Swap rates by a boot-strapping procedure
%--------------------------------------------------------------------------
% Array of discount Factors at each payment date on fixed side 
discountFactors = []; 

tic; 
for i = 1:length(swapRates), 
   
% Use all the pre-calculated forward rate data in boot-strapping procedure: 
if ~isempty(accrualDateNumbers),
    startIndex = max(find(datenum(swapAccrualDates) <= accrualDateNumbers(end)));
    fprintf('\n S_j: %s T_l: %s\n',swapAccrualDates{startIndex}, datestr(accrualDateNumbers(end)));
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
    fprintf(' %s to %s: index: %d\n',datestr(accrualDateNumbers(end-1)),swapAccrualDates{j+1},datenum(swapAccrualDates(j+1)) - accrualDateNumbers(end-1));     
end

forwardRates(end+1) = Bisection(M, 0.0, 1.0, 0.000000001, F);

% interpolate of the forward curve assuming *constant* forward rates  
if length(forwardRates) > 1, 
    forwardCurve((accrualDateNumbers(end-1)+1 - datenum(startDate)):(accrualDateNumbers(end) - datenum(startDate))) = forwardRates(end);
else
    forwardCurve(1:(accrualDateNumbers(end) - datenum(startDate))) = forwardRates(end);
end

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

% Plot forward rates derived from Deposit rates and futures and Swaps 
if 0
    figure('Position',[10 300 600 400])
    box on;
    hold on;
    t = accrualDateNumbers - datenum(startDate);
    indices = 1:length(deposits); 
    plot(t(indices)/365,forwardRates(indices),'rx')
    indices = (indices(end)+1):(indices(end) + length(futures)); 
    plot(t(indices)/365,forwardRates(indices),'r*')
    indices = (indices(end)+1):(indices(end) + length(swapRates)); 
    plot(t(indices)/365,forwardRates(indices),'ro')
    plot((1:t(end))/365,forwardCurve,'b-')
    
    % Compute discount curve from forward curve
    discountCurve = zeros(length(forwardCurve),1); 
    for i = 1:length(forwardCurve)
        discountCurve(i) = exp(-sum((1/365)*forwardCurve(1:i))); 
    end

    % compute the yield Curve or YTM from the discount Curve  
    yieldCurve = -log(discountCurve)' ./ ((1:length(discountCurve))/365); 
       
    % 2-D line plots with y-axes on both left and right side
    [AX, H1, H2] = plotyy((1:t(end))/365,yieldCurve,(1:t(end))/365,discountCurve,'plot','bar');
        
    set(H1,'LineStyle','-','LineWidth',1.0,'Color',[0.25 0.25 0.25]);
    set(get(AX(1),'Ylabel'),'String','Forward rate (instantaneous)','FontSize',12);
    set(AX(1),'YColor','k');
    xlim(AX(1),[0 swapRates(end,1)])
    ylim(AX(1),[0 0.08])
    set(H2,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
    set(get(AX(2),'Ylabel'),'String','Discount Curve','FontSize',12);
    set(AX(2),'YColor','k');
    xlim(AX(2),[0 swapRates(end,1)])
    
    % 
    % when permutating the handles need to change the color of the front axes to "none"
    %
    fc = get(gcf, 'Children');
    set(AX(1), 'Color', 'none');
    set(AX(2), 'Color', 'w');
    set(gcf, 'Children', flipud(fc));
        
    xlabel('Benchmark maturities (red)','FontSize',12,'Color','k'); 
    title('The instantaneous forward rate f(0,t)','FontSize',14,'Color','k');   
end  

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
    
    indices = 1:length(deposits); 
    t = accrualDateNumbers - datenum(startDate);
    % plot(t(indices)/365,forwardRates(indices),'rx')
    if ~isempty(indices), 
        plot(t(indices)/365,forwardRates(indices),'kx','LineWidth',0.5)
        indices = (indices(end)+1):(indices(end) + length(futures)); 
    else 
        indices = 1:length(futures); 
    end
    
    if ~isempty(indices),    
        plot(t(indices)/365,forwardRates(indices),'kx','LineWidth',0.5)
        indices = (indices(end)+1):(indices(end) + length(swapRates)); 
    else
        indices = 1:length(swapRates);
    end
    
    % 
    % plot(t(indices)/365,forwardRates(indices),'ro')
    plot(t(indices)/365,forwardRates(indices),'ko','LineWidth',0.75)
    % plot((1:t(end))/365,forwardCurve,'b-')
    %%% plot((1:t(end))/365,forwardCurve,'k-','LineWidth',1.0)
    maturity = (1:t(end))/365;
    stairs(maturity,forwardCurve,'k-')
    
    plot(maturity,yieldCurve,'k-','LineWidth',0.75)
    
    %ylabel('Forward Rate Curve','FontSize',14,'Color','k');
    ylabel('Rate','FontSize',14,'Color','k');
    
    % extract year from file name (for this to work we must ensure file-name syntax is consistent)
    year = file_name(16:19);
    
    if strcmp(year,'2011'), 
        axis([0 30.02 0 0.06]) 
    elseif strcmp(year,'2007'),
        axis([0 30.02 0.04 0.06])  
    end
        
    box on;
    hold on; 
    xlabel('Benchmark maturities','FontSize',14,'Color','k');   
    %text(15,0.05,'\it forward-rate curve','FontSize',14')  
    %title('The Instantaneous Forward Rate f(t,T)','FontSize',16,'Color','k');   
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
     plot((datenum(swapAccrualDates(2*swapRates(i,1) + 1)) - datenum(startDate)),discountFactors(i),'ko')
     hold on; 
end
plot(discountCurve,'b-');

%--------------------------------------------------------------------------
%  Data output for Vasicek short-rate calibration exercise 
%--------------------------------------------------------------------------
if IS_VASICEK
    % approximation to the instantaneous short rate r(0):  
    fprintf('O/N rate: 4.2f\n',deposits{1,2}); 
    
    shortRate = 0.01*deposits{1,2}; 
    
    t = accrualDateNumbers - datenum(startDate);
    % Output discount factors for each year starting at 'startDate'
    for i = 1:10,     
        index = datenum(datevec(startDate) + [i 0 0 0 0 0]); 
        fprintf(' %s & %d & %2.4f \n', datestr(index),i,discountCurve(index -(datenum(startDate) + t(1)) + 1));

        vasicek(i) = discountCurve(index -(datenum(startDate) + t(1)) + 1);
    end

    save('vasicek_data.mat', 'vasicek', 'shortRate');    
    %--------------------------------------------------------------------------
    % LMFSOLVE  Solve a Set of Nonlinear Equations via Least-Squares.
    % Levenberg-Maquardt used to minimise a sum of squared residuals. 
    %
    % beta = [    ] for an intial guess of [0.1 0.1 0.1]
    %
    %--------------------------------------------------------------------------

    %path(path,'C:\MATLAB_SV701\work\LMFsolve[1]')
    path(path,'C:\Users\ania\Documents\SIR\code\LMFsolve[1]');
    theta = 0.001;
    alpha = 0.05;
    sigma = 0.015; 
     
    %theta: 0.02200
    %alpha: 0.46507
    %sigma: 0.01229
    
    beta  = [theta alpha sigma]';
    beta  = [0.05*rand rand 1.5*rand]';  
    % [beta,ssq,cnt] = LMFsolve('vasicek_ex1r',beta);
    
    opts = [1 1e-8 1e-8 1000 0.25*1e-5];
    tic;
    [beta,info,perf] = SMarquardt('vasicek_ex1r',[],beta,opts);
    % [theta,info,perf] = marquardt('correlation_ex3r',[],beta,opts(1:4));
    t = toc; 

    fprintf(' No of iterations:         % d\n', info(5)); 
    fprintf(' Sum of squared residuals: % 1.4f\n',info(1));
    fprintf(' Inf norm of gradient:     % 2.2e\n',info(2)); 
    fprintf(' Norm of dx:               % 2.2e \n',info(3));
    fprintf(' Time taken (seconds):     % 2.2f \n\n',t); 

    % beta = [0.022 0.46507 0.01229]';
    fprintf(' theta: %5.5f\n',beta(1,end)); 
    fprintf(' alpha: %5.5f\n',beta(2,end));     
    fprintf(' sigma: %5.5f\n',beta(3,end));
    fprintf(' theta/alpha: %5.5f\n',beta(1,end)/beta(2,end)); 
    
    % Constrain value to be at least 1% 
    if beta(3,end) < 0.01, 
        beta(3,end) = 0.01; 
        fprintf(' Apply Constraint that value of beta is to be at least 1 per cent \n'); 
    end
    
    modelBondPrice = vasicek - vasicek_ex1r(beta(:,end),[])';
    Error          = vasicek - modelBondPrice ./ vasicek; 
    
    % Output data in the form of a table with %-errors 
    for i = 1:length(modelBondPrice), 
        fprintf(' %d & %1.4f & %1.4f & %1.2f \\\\ \n',i,vasicek(i), modelBondPrice(i), Error(i));    
    end
    
    if 1
    figure;
    plot(1:10,modelBondPrice','k:o'); 
    hold on; 
    plot(1:10,vasicek,'k-*'); 
    
    ylabel('Bond price','FontSize',14,'Color','k');
    box on;
    xlabel('Maturities','FontSize',14,'Color','k');     
    %title('Model and Market Bond Prices','FontSize',16,'Color','k');  
    
    end
end
    
return;