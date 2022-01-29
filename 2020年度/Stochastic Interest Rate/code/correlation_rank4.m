%--------------------------------------------------------------------------
% Examples: 
% The most general methodology to create a valid correlation
% matrix for risk management and option pricing purposes 
% rank-4 matrix 
%--------------------------------------------------------------------------
clear all;
%--------------------------------------------------------------------------
% SMarquardt: Solve a Set of Nonlinear Equations via Least-Squares.
% Levenberg-Maquardt used to minimise a sum of squared residuals. 
%--------------------------------------------------------------------------
path(path,strcat(pwd,'\LMFsolve[1]'));

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
 



