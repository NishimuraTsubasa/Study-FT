%--------------------------------------------------------------------------
% Examples: 
% The most general methodology to create a valid correlation
% matrix for risk management and option pricing purposes 
% rank-4 matrix 
%--------------------------------------------------------------------------
clear all;
if 0 
if 0 
LongCorr = 0.5;
beta = 0.05;
N =10; 

for i = 1:N,
  for j = 1:N, 
    rho(i,j) = LongCorr + (1-LongCorr)*exp(-beta*abs(i-j));
  end
end
end

N = 20; 
load doust

[V,D] = eig(rho);

% Order the eigenvalues from highest to lowest 
[E,I] = sort(diag(D),1,'ascend');

% Create associated matrix of eigenvectors (orthonormal set: U*U' = I)
U = V(:,I);

% flag eigenvalues less than some pre-defined value  
% I = find(E<0.04);

% 4-factor model -> still retains approximate shape of sample correlation matrix
%I = find(E<0.15); 
I = find(E<0.2);

I = I(end)+1;

U = U(:,I:length(E)); 
E = E(I:length(E)); 

% Matrix 'B'; pseudo-sqaure root of the correlation matrix
B = U*diag(sqrt(E));

% The diagonal scaling matrix T wrt the eigensystem: 
% normalising the row vectors of B to unit length ...
% necessary if focusing only on the first few principal components
% or the eigenvectors corresponding to positive eigenvalues 
for i = 1:length(U),
  T(i) = 1/(sum((U(i,:).^2).*transpose(E)));
end
T = diag(T);

% Reduced rank correlation matrix
C = sqrt(T)*B*(sqrt(T)*B)';

% or normalise by treating B*B' as a covariance matrix 
C = B*B';
C = C ./(sqrt(diag(C) * diag(C)'));

% Note that: B'B = diag(E) 
%            (sqrt(T)*B)'*(sqrt(T)*B) = {sum(diag(T) .* U(:,i) .* U(:,j))*sqrt(E(i)*E(j))} 
%                                       which for a given T does not form a diagonal matrix  
% Determine reduced rank matrix using angles parameterisation and 
% then orthogonalise new matrix                                        

%--------------------------------------------------------------------------
% By inverting the relevant transformation we can find the angles theta: 
% Invert 'non-normalised B' to get 
%  theta = [1.6695 1.7239 1.2837; ... 
%           1.5914 1.6672 1.3068; ...  
%                               ];
%--------------------------------------------------------------------------
 
for i = 1:length(B),
  theta(i,1) = acos(B(i,1)); 
  theta(i,2) = acos(B(i,2)/sin(theta(i,1)));
  theta(i,3) = acos(B(i,3)/(sin(theta(i,1))*sin(theta(i,2))));  
  theta(i,4) = asin(B(i,4)/(sin(theta(i,1))*sin(theta(i,2)))); 
% should equal theta(i,3), not the case here - only for normalised B will this be 'true' 
% ('true' -> can only get first few entries to match) 
end

% Test to see if we recover correct reduced rank correlation matrix
%for i = 1:length(B),
%  B(i,:) = [cos(theta(i,1)), cos(theta(i,2))*sin(theta(i,1)), cos(theta(i,3))*sin(theta(i,1))*sin(theta(i,2)), ... 
%                                                                   sin(theta(i,1))*sin(theta(i,2))*sin(theta(i,3))];
%end 
%C = B*B'; % --- must invert normalised B to get orginal correlation matrix

%--------------------------------------------------------------------------
% Theta via numerical optimisation 
%--------------------------------------------------------------------------
if 0 
theta = [1.6844  1.7328 1.2775;
         1.6088  1.6828 1.2965;
         1.4688  1.581  1.3444;
         1.4435  1.4708 1.4267;
         1.5051  1.3957 1.5203;
         1.6365  1.3957 1.6213;
         1.6981  1.4708 1.7149;
         1.6728  1.581  1.7972;
         1.5328  1.6828 1.8451;
         1.4571  1.7328 1.864];
end

% resulting optimal rank-4 correlation matrix
for i = 1:length(B),
  B(i,:) = [cos(theta(i,1)), cos(theta(i,2))*sin(theta(i,1)), cos(theta(i,3))*sin(theta(i,1))*sin(theta(i,2)), ... 
                                                              sin(theta(i,1))*sin(theta(i,2))*sin(theta(i,3))];
end 
CTheta = B*B'; 

% note sigmoid-like shape due to reduced rank approximations
%for i=1:length(rho),
%  hold on;  
%  plot(rho(:,i),'-+');
%  plot(C(:,i),':d');   
%  plot(CTheta(:,i),'--o');
%  pause; 
%  close;
%end
end
%--------------------------------------------------------------------------
% SMarquardt: Solve a Set of Nonlinear Equations via Least-Squares.
% Levenberg-Maquardt used to minimise a sum of squared residuals. 
%--------------------------------------------------------------------------
% path(path,'C:\Program Files\Microsoft Visual Studio 9.0\VC\bin')

%path(path,'C:\MATLAB_SV701\work\LMFsolve[1]')
path(path,'C:\Users\ania\Documents\SIR\code\LMFsolve[1]');

N = 20;
x0 = rand([1,3*N]); 
opts = [1 1e-11 1e-10 1000 0.25*1e-5];
tic;
[theta,info,perf] = SMarquardt('correlation_ex3rb',[],x0,opts);
%[theta,info,perf] = marquardt('correlation_ex3r',[],x0,opts(1:4));
t = toc; 
fprintf(' No of iterations:         % d\n', info(5)); 
fprintf(' Sum of squared residuals: % 1.4f\n',info(1));
fprintf(' Inf norm of gradient:     % 2.2e\n',info(2)); 
fprintf(' Norm of dx:               % 2.2e \n',info(3));
fprintf(' Time taken (seconds):     % 2.2f \n\n',t); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When using the analytical Jacobian with the above opts(1:4) we get warning: 
%
% Warning: Matrix is close to singular or badly scaled.
%          Results may be inaccurate. RCOND = 1.866770e-016.
%
% For the opts chosen the damping parameter mu can become very small ~ 1e-15
% and the matrix becomes singular. 
% Changed marquardt to ensure that mu is not less than 1e-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Min damping parameter:    % 2.2e \n',min(perf(3,:))); 

theta = mod(theta(:,end),2*pi); 
%fprintf(' Theta(0): % 1.4f Theta: % 1.4f\n', [x0', theta]');  
%fprintf(' No of iterations: %d\n Sum of squared residuals: % 1.4f\n\n',cnt,ssq); 

theta = reshape(theta,3,N)';
for i = 1:N,
  B(i,:) = [cos(theta(i,1)), cos(theta(i,2))*sin(theta(i,1)), cos(theta(i,3))*sin(theta(i,1))*sin(theta(i,2)), ... 
                                                              sin(theta(i,1))*sin(theta(i,2))*sin(theta(i,3))];
end 
C = B*B'; % as given on page 257? 
%fprintf(' Norm of (C-CTheta):       % 2.2e \n\n',norm(C-CTheta));         
          
% The columns of B are not orthogonal -> diagonalise rank-4 correlation matrix 
%[B, U, E, V, D] = pseudosquareroot(C, 10^-8);

%disp([' Matrix \hat{B}: ']); 
%fprintf('% 1.5f % 1.5f \n',B')
%disp([' Matrix \hat{C}: ']); 
%fprintf('% 1.4f % 1.4f % 1.4f % 1.4f % 1.4f % 1.4f % 1.4f % 1.4f % 1.4f % 1.4f\n',C)
%disp([' '])

%disp([' Test to see if columns of \hat{B} are orthogonal: ', num2str(sum(B(:,1) .* B(:,2)))]);

colormap gray  
mesh(C)
axis([1 20 1 20 0.3 1])
zlabel(['Correlation'],'FontSize',14,'Color','k'); 
xlabel(['Forward-rate expiry'],'FontSize',14,'Color','k');
ylabel(['Forward-rate expiry'],'FontSize',14,'Color','k');
%title(['Correlation Matrix (4-factor model)'],'FontSize',14,'Color','k');
view(-30,30);

% Make up one column of the colormap. Now the colormap goes from 0 at index 1
% corresponding to gray level 24, to 1 at index 256 corresponding to gray level 256.
cmap1 = [linspace(0, 0.3, 256-0)];
% Make up a gray colormap by copying this to all 3 columns.
cmap3 = [cmap1; cmap1; cmap1]';
% Apply the colormap to the image.
colormap(cmap3);
 



