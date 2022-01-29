%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daragh McInerney, 2012
%
% Probability density of CEV and displaced diffusion process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

sigma     = 0.1;
beta      = 0.5;   
T         = 1.0; 
N         = 50;
timeStep  = T/N;

% use these settings for F(0) = 1% 
F0       = 0.01;
stepSize = 0.0005;
% use these settings for F(0) = 5% 
F0       = 0.05;
stepSize = 0.002;

%-------------------------------------------------------------------------- 
% Distribution of forward rate is noncentral chi-squared in CEV model. 
% Simulate the distribution of F at time T: 
%--------------------------------------------------------------------------
M      = 100000; % number of simulation trials 
F      = zeros([M,N]); 
F(:,1) = F0; 
StDev  = sigma*sqrt(timeStep);
count  = 0; 

for i = 1:M, 
  random = randn(N,1);
  for j=2:N, 
    F(i,j) = F(i,j-1) + StDev*(F(i,j-1)^beta)*random(j); 
    if (F(i,j) < 0)                                  % If F becomes negative -> process no longer well-defined 
      F(i,j:N) = 0.0;                                % ... Set remaining values to zero    
      count = count + 1;
      break;
    end    
  end
end
fprintf('\n Percentage of CEV processes that hit absorbing boundary at zero: %0.5f\n',count/M); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact Monte-Carlo of CEV process
%
% I don't know a direct reference for that, but the stock price distribution is known in terms
% of the noncentral chisquared distribution (see Schroder's paper on computing the CEV model).
% Then, you only need to know how to draw from that, which should be discussed in Glasserman or
% Brodie & Kaya's 'exact simulation' paper, or elsewhere.
%
% p.s. Thinking harder, the stock price distribution for S(T) should be partitioned into two event sets:
% i) absorption at S=0 has not occurred prior to T; ii) absorption has occurred.
% Then, you get i) by drawing from the nc chisquared. You get ii) by drawing from the exact
% absorption distribution, which is a gamma distribution. Alternatively, if you don't want to simulate
% ii), just use the exact formula for the absorption probability: P_abs(S0,T) and 
% your simulation (schematically) is:
%
% S(T) = P_abs(S0,T) DiracDelta(S=0) + (1 - P_abs(S0,T)) (Draw from nc chisquared)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eqn (1) in Schroder gives the density for S(T), valid for all S(T) > 0.
% This is all you need to simulate the call value. 
%
%p.s.
%------------------------------------------------------------------------------------------
%Something like Mathematica would help, but 
%if whatever platform you use has the Bessel function of eqn(1), try coding this (again schematic):
%
%Code F(ST) = int_0^ST pdf(x) dx, where pdf(x) = eqn(1) from Schroder
%Draw ST using the Inverse Method (given a continuous uniform variable U in [0,1] and an invertible 
%cumulative distribution function F, the random variable X = F^{?1}(U) has distribution F).
%Now, you should first confirm that F(x) is a defective distribution: F(infinity) = 1 - P_abs(S0,T) < 1.
%That means that when you make a random draw U such that F(infinity) < U <= 1, then there will be no root.
%To do this 'bad draw' check, I suggest working out F(infinity) analytically -- it's not too hard.
%Then, just interpret the bad draws as {S(T) = 0}, hence C(T) = 0, so continue() to the next draw (without a root search) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here we assume that max(F(:,N)) <= 1 
x = 0:stepSize:1; 
% find histrogram to show distribution of data values
y = hist(F(:,N),x)/(M*stepSize); 

%--------------------------------------------------------------------------
% For initial values of the forward rate close to zero (0.005) and a high 
% enough value of sigma the pdf is given by a delta function at F = 0 
%--------------------------------------------------------------------------
box on;
hold on; 
bar(x,y,1)
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12')

% Change the color of the graph so that the bins are white and the edges 
% of the bins are black.    
h = findobj(gca,'Type','patch');
%set(h,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]);   % 'Color',[0 0 1]
set(h,'FaceColor',[1.0 1.0 1.0],'EdgeColor',[0.0 0.0 0.0]);

[X, I] = sort(F(:,N)); 
axis([-stepSize mean(X(end-(M/500):end)) 0 floor(max(y)+0.1*max(y))]);

%--------------------------------------------------------------------------
% Density of CEV process
% To calculate The density function of F(T) conditional on the current forward F(0)
% we use Equation (1) from M. Schroders paper: 
% “Computing the Constant Elasticity of Variance Option Pricing Formula,” Journal of
% Finance, v. 44 (March 1989), 211-219.
%--------------------------------------------------------------------------
% nc chisquared distribution parameters   
q     = 0.0; 
r     = 1e-6; % function fails for r = 0 so need to set r to some small value
beta2 = 2*beta; 

k = (2*(r-q))/((sigma^2)*(2-beta2)*(exp((r-q)*(2-beta2)*T)-1));
x = k * ((F(1))^(2-beta2))*exp((r-q)*(2-beta2)*T);

z = 0.0:stepSize:1;
w = k * z.^(2-beta2);
f = (2-beta2)*k^(1/(2-beta2))*(x*w.^(1-2*beta2)).^(1/(4-2*beta2)).*exp(-x-w).*besseli(1/(2-beta2),2*sqrt(x.*w));

% If second argument in bessel function is large ~ 1000 then value is Inf 
hold on 
box on; 
plot(z,f,'k','LineWidth',1); 
%title('Density of Forward Rate at Time T = 1','FontSize',16,'Color','k')
ylabel('Probability density function','FontSize',16,'Color','k')
xlabel('Forward rate','FontSize',16,'Color','k') 


%--------------------------------------------------------------------------
% Model the forward rate process as a Displaced-Diffusion 
% Simulate the distribution of F at time T: 
%--------------------------------------------------------------------------
if 1
% Parameters for displaced diffusion process in terms of CEV process 
sigmaHat = sigma*beta*(F0^(beta-1));
eta      = (1- beta)*F0/beta;

M = 100000; % number of simulation trials  
F = zeros([M,N]); 
F(:,1) = F0; 
for i = 1:M, 
  random = randn(N,1);
  for j=2:N, 
      % Note: F > -eta as F + eta is lognormal and therefore greater than 0
      F(i,j) = F(i,j-1) + sigmaHat*(F(i,j-1) + eta)*sqrt(timeStep)*random(j);  
  end
end

%--------------------------------------------------------------------------
% Compare with log-normal distribution ... Calculate mean and standard 
% deviation of the variable's natural logarithm
%--------------------------------------------------------------------------
x = 0:stepSize:1;
% mean of log F(T) (exact) 
mu = log(F0 + eta) - 0.5*(sigmaHat^2)*T; 
% standard deviation of log F(T) (exact)
StDev = sigmaHat*sqrt(T);                    

% aggrement between the vaules calculated using simulation and the exact formula  
fprintf('mean of log (F(T) + eta):  % 1.5f (simulation) % 1.5f (exact)\n',mean(log(F(:,N) + eta)),mu);
fprintf('StDev of log (F(T) + eta): % 1.5f (simulation) % 1.5f (exact)\n',std(log(F(:,N) + eta)),StDev);

logNormalDensity = (1./(sqrt(2*pi)*StDev.*(x))).*exp(-((log(x) - mu).^2)./(2*StDev^2)); 
figure;

set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12')



hold on;
box on; 
%plot(zeros(1,2),0:100:100,'k:','LineWidth',2); 
plot(x-eta,logNormalDensity,'k-');
plot(x,f,'k--'); 
axis([-0.01 0.15 0 floor(max(f)+0.1*max(f))]);
%title('Density of Forward Rate at Time T = 1','FontSize',16,'Color','k')
ylabel('Probability density functions ','FontSize',16,'Color','k')
xlabel('Forward rate','FontSize',16,'Color','k') 
end 


