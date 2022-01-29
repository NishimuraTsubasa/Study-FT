%--------------------------------------------------------------------------
% Ploting Doust Correlation structure for forward rates with convexity
%--------------------------------------------------------------------------

clear all;
N          = 20;       % number of forward rates 
delta_t    = 1;        % accrual period 
time_steps = [1:delta_t:N];

%--------------------------------------------------------------------------
path(path,strcat(pwd,'\LMFsolve[1]'));
x = rand([1,2]);
%x = rand([1,5]);
[Beta,ssq,cnt]=LMFsolve('doust_correlation_empiricalr',x,'MaxIter',5000)

fprintf(' Beta(0): % 1.4f Beta: % 1.4f\n', [x', Beta]');  
fprintf(' No of iterations: %d\n Sum of squared residuals: % 1.4f\n\n',cnt,ssq); 

Beta   = Beta(1) * (1./((1:(N-1)).^Beta(2)));
%Beta   = Beta(1) + Beta(2) * (1./((1:(N-1)).^1)) + Beta(3) * (1./((1:(N-1)).^2)) ... 
%       + Beta(4) * (1./((1:(N-1)).^3)) + + Beta(5) * (1./((1:(N-1)).^4)); 
a      = exp(-Beta*delta_t);

sum = 1; 
for i = 1:(length(Beta)+1),      
  for j = 1:(length(Beta)+1),
     
     if (j < i), 
       for k = j:(i-1),
         sum = sum * a(k);
       end 
       rho(i,j) = sum; 
    
       rho(j,i) = rho(i,j);
       sum = 1;
     end   
  end
  rho(i,i) = 1; 
end

if 1 
    colormap gray   
    mesh(rho)
    hold on;
    axis([1 20 1 20 0.3 1])
    zlabel(['Correlation'],'FontSize',14,'Color','k'); 
    xlabel(['Forward-rate expiry'],'FontSize',14,'Color','k');
    ylabel(['Forward-rate expiry'],'FontSize',14,'Color','k');
    %title(['Correlation Matrix'],'FontSize',14,'Color','k');
    view(-40,30);

    % Make up one column of the colormap. Now the colormap goes from 0 at index 1
    % corresponding to gray level 24, to 1 at index 256 corresponding to gray level 256.
    cmap1 = [linspace(0, 0.3, 256-0)];
    % Make up a gray colormap by copying this to all 3 columns.
    cmap3 = [cmap1; cmap1; cmap1]';
    % Apply the colormap to the image.
    colormap(cmap3);
end

save('doust_Mar2013','rho'); 