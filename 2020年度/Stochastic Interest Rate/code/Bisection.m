function [x] = Bisection(Target, Low, High, Tolerance, F)

% Check to see if Function is monotonically increasing or decreasing 
if (F(High) > F(Low)),  
   sign = 1; 
else
   sign = -1;
end 

x = 0.5*(Low + High);
y = F(x);

count = 1; 

while (abs(y-Target) > Tolerance)
   
    if (sign*(y - Target) < 0)
        Low = x;
    else
        High = x;
    end

    x = 0.5*(Low+High);
    y = F(x);
    
    count = count + 1; 
    
    %if count > 100000,
    %    fprintf('Exited with accuracy of %1.8f\n', abs(y-Target)); 
    %    break
    %end
end

return;