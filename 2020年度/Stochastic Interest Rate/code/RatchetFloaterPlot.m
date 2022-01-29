clear all;
format long;

% load ratchet cash-flows data from a csv file
file_name = 'RatchetFloater_Table.csv';

fid = fopen(file_name,'r');

% Set headerLines to 1 to skip the first line in the file 
C = textscan(fid,'%f%f%f%f%f%f%*[^\n]','delimiter',',','headerLines',1);

TimeToExpiry = C{1}

for j = 1:length(C)-2,
    for i = 1:length(C{1}) 
        cashFlows(i,j) = C{j+2}(i)
    end 
end
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12')

plot(TimeToExpiry,cashFlows,'kx-')
ylabel('Discounted cash flows','FontSize',16,'Color','k');
xlabel('Cash-flow times','FontSize',16,'Color','k');  
axis([0.5 3.0 -200 1800])

text(2.5,1550,'$\alpha = 0.25$','Interpreter','latex','FontSize',14)
text(2.5,1000,'$\alpha = 0.1$','Interpreter','latex','FontSize',14)
text(2.5,400,'$\alpha = 0.05$','Interpreter','latex','FontSize',14)
text(2.5,-50,'$\alpha = 0.03$','Interpreter','latex','FontSize',14)