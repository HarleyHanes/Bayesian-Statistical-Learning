%-------------------------------------------------------------------
% Harley Hanes                                               8/29/19
% Bayesian Statistical Learning                           Homework 1
%-------------------------------------------------------------------
clear; close all;
%Data Input and Basis Function Declaration
x_data=[0; 1; 2; 3; 4];
y_data=[1.09; .5; -.94; -3.57; -7.02];
% Section 1
rho1 = @(x) 0*x+1;
% Section 2
rho2 = @(x)[0*x+1 x];
% Section 3
rho3 = @(x)[0*x+1 x x.^2];
% Section 4
rho4 = @(x)[0*x+1 x x.^2 x.^3];
BasisFunctions={rho1 rho2 rho3 rho4};

% Coeffecient Calculation
Beta=zeros(4);
for i=1:4
    rho=BasisFunctions{i};
    X=rho(x_data);
    Beta(1:i,i)=(X'*X)\(X'*y_data);
end

%Plotting
plot(x_data,y_data,'r*','MarkerSize',10)
hold on 

f = @(Coef, x)Coef(1) + Coef(2)*x +Coef(3)*x.^2 +Coef(4)*x.^3;
x=-1:.01:5;
for i=1:4
    line=f(Beta(:,i),x);
    plot(x,line,'LineWidth',1.5)
end
labels{1}='Data Points';
labels{2}=sprintf('\\rho (x)= %.2f',Beta(1,1));
labels{3}=sprintf('\\rho (x)= %.2f + %.2fx',Beta(1:2,2));
labels{4}=sprintf('\\rho (x)= %.2f + %.2fx + %.2fx^2', Beta(1:3,3));
labels{5}=sprintf('\\rho (x)= %.2f + %.2fx + %.2fx^2 + %.2fx^4', Beta(1:4,4));
legend(labels);

xlabel('x')
ylabel('y')

%% lsqnonlin and fminsearch practice
clear;
x=[0; 1; 2; 3; 4];
y_data=[1.09; .5; -.94; -3.57; -7.02];
f = @(Coef)Coef(1) + Coef(2)*x +Coef(3)*x.^2 +Coef(4)*x.^3-y_data;
init=[1 0 -.6 .1];

lsqnonlin(f,init)

lsq= @(Coef)sum((f(Coef)-y_data).^2);
fminsearch(lsq,init)


