%===================================================================
% Harley Hanes
% Functions for Solving Linear Systems                      
% Bayesian Statistical Learning                           
%-------------------------------------------------------------------
% This is the main file for the code set to solve linear systems using 
% linear and nonlinear solvers. Code is broken into following
% sections.
% 1) Master control- Solver, data, and basis function selection
% 2) Data and basis function entry
% 3) Liner solver
% 4) Lsqnonlin
% 5) fminsearch
% 6) Plotting
%-------------------------------------------------------------------
% Version Notes
% 1)Haven't figured out how to adjust input vector lengths based on # of
%   basis functions. Is creating issue for fminsearch
% 2)fminsearch seems to need a decrease in tolerance to work for the class
%   data
% 3)Plotting currently only works for Linear Solve and legend doesn't work
%===================================================================
clear; close all;clc;
%% Master Control
DataType='Class Data'; %Options: 'Class Data', 'Normal Scatter
                          %         'Straight Line'
BasisFunc='3rd Order Poly'; %Opions '0, 1st, 2nd, 3rd Order Poly'
Solver='all'; %Options: Linear, lsqnonlin, fminsearch, all
PlotBool='yes';


%% Data and basis function entry
switch DataType
    case 'Class Data'
        xData=[0; 1; 2; 3; 4];
        yData=[1.09; .5; -.94; -3.57; -7.02];
    case 'Normal Scatter'
        xData=(1:.1:10)';
        yData=(randn(1,length(xData))*.1)';
    case 'Straight Line'
        xData=(1:.1:10)';
        yData=ones(1,length(xData))';
    otherwise 
        fprintf('Error!! DataType not recognized')
        keyboard
end

switch BasisFunc
    case '0 Order Poly'
        rho=@(x) x.^0;
        fModel = @(Coef)Coef(1);
    case '1st Order Poly'
        rho=@(x)[x.^0 x.^1];
        fModel = @(Coef)Coef(1) + Coef(2)*xData;
    case '2nd Order Poly'
        rho=@(x)[x.^0 x.^1 x.^2];
        fModel = @(Coef)Coef(1) + Coef(2)*xData +Coef(3)*xData.^2;
    case '3rd Order Poly'
        rho=@(x)[x.^0 x.^1 x.^2 x.^3];
        fModel = @(Coef)Coef(1) + Coef(2)*xData +Coef(3)*xData.^2 +Coef(4)*xData.^3;
     otherwise
        fprintf('Error!! BasisFunc not recognized')
        keyboard
end
%% Linear Solver
X=rho(xData);
LinBeta=((X'*X)\(X'*yData))';

%% Lsqnonlin

lsqfunc= @(Coef)fModel(Coef)-yData;
init=zeros(1,4);
lsqBeta=lsqnonlin(lsqfunc,init);

%% fminsearch
fR=@(Coef)fModel(Coef)-yData;
fminfunc= @(Coef)fR(Coef)'*fR(Coef);
options = optimset('TolFun',1e-8,'TolX',1e-8);

init=zeros(1,4);
fminBeta=fminsearch(fminfunc,init,options);
%% Display and Plotting
switch Solver
    case 'all'
        DispLinBeta=sprintf('%.2f,',LinBeta);
        DisplsqBeta=sprintf('%.2f,',lsqBeta);
        DispfminBeta=sprintf('%.2f,',fminBeta);
        DispMessage=sprintf('Model Coefficients for each Solver\nLinear: %s\nlsqnonlin: %s\nfminsearch: %s',DispLinBeta,DisplsqBeta,DispfminBeta);
        disp(DispMessage)
    otherwise
        fprintf('Error!! Solver not recognized')
        keyboard
end
if strcmpi(PlotBool,'yes') || strcmpi(PlotBool,'no')
    tPlot=linspace(min(xData)-(min(xData)+max(xData))/10,...
        max(xData)+(min(xData)+max(xData))/10);
    yPlot=zeros(1,length(tPlot));
    if size(LinBeta)== size(rho(xData(1)))
        Beta=LinBeta';
    else
        Beta=LinBeta;
    end
    for i=1:length(tPlot)
        yPlot(i)=rho(tPlot(i))*Beta;
    end
    plot(xData,yData,'r*','MarkerSize',10)
    hold on 
    plot(tPlot,yPlot)
    PlotLegend={'Data',sprintf('%s\nBeta=%.2f',BasisFunc,Beta')};
    legend(PlotLegend)
    xlabel('x')
    ylabel('y')
    PlotTitle=sprintf('%s Function Fit for %s',BasisFunc,DataType);
    title(PlotTitle)
end



