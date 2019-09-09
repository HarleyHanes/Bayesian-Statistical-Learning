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
% 1)fminsearch is not reaching its optimal value for Runge function (giving
% a sum of squared residuals of approximately .02
% 2)ODE function isnt minimizing because each solve iteration has different
% numbers of data points so I can't line up data
%===================================================================
clear; close all;clc;
%% Master Control
DataType='Harmonic ODE'; %Options: 'Class Data', 'Normal Scatter
                          %         'Straight Line', 'Runge Function'
BasisFunc='Harmonic ODE'; %Opions '0, 1st, 2nd, 3rd Order Poly'
                            %        'Mixed Exponential'
                            %        'Self Enter
Solver='fminsearch'; %Options: Linear, lsqnonlin, fminsearch, all
PlotBool='yes';


%% Data and basis function entry
%Note: Data vectors need to be column vectors for fminsearch and linear
%Solver
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
    case 'Runge Function'
        NumPoints=20;
        xmin=-1;
        xmax=1;
        a=1;
        b=1;
        c=25;
        TrueFunc=@(x) a./(b+c.*x.^2);
        xData=linspace(xmin,xmax,NumPoints)';
        yData=TrueFunc(xData);
    case 'Harmonic ODE'
        NumPoints=20;
        tspan=[0,10];
        Coef=[0 1/2 1 2];
        numCoef=length(Coef);
        [xData,yData]=HarmonicODE(Coef,tspan,NumPoints);
    otherwise 
        fprintf('Error!! DataType not recognized')
        keyboard
end

switch BasisFunc
    case '0 Order Poly'
        rho=@(x) x.^0;
        fModel = @(Coef)Coef(1);
        numCoef=1;
    case '1st Order Poly'
        rho=@(x)[x.^0 x.^1];
        fModel = @(Coef)Coef(1) + Coef(2)*x;
        numCoef=2;
    case '2nd Order Poly'
        rho=@(x)[x.^0 x.^1 x.^2];
        fModel = @(Coef)Coef(1) + Coef(2)*x +Coef(3)*x.^2;
        numCoef=3;
    case '3rd Order Poly'
        rho=@(x)[x.^0 x.^1 x.^2 x.^3];
        fModel = @(Coef,x)Coef(1) + Coef(2)*x +Coef(3)*x.^2 +Coef(4)*x.^3;
        numCoef=4;
    case 'Mixed Exponential'
        fModel = @(Coef,x)Coef(1) + Coef(2) * exp(Coef(3).*x + Coef(4)*x.^2);
        numCoef=4;
    case 'Harmonic ODE'
        tspan=[0,10];
        fModel=@(Coef)HarmonicODE(Coef,tspan,NumPoints);
    case 'Self Enter'
    otherwise
        fprintf('Error!! BasisFunc not recognized')
        keyboard
end
%% Linear Solver
if strcmpi(Solver,'Linear') || strcmpi(Solver,'All')
    X=rho(xData);
    linBeta=((X'*X)\(X'*yData))';
    if size(linBeta)== size(rho(xData(1)))
        linBeta=linBeta';
    end
    if strcmpi(Solver,'linear')
        Beta=linBeta;
    end
end

%% Lsqnonlin
if strcmpi(Solver,'Lsqnonlin') || strcmpi(Solver,'All')
    if strcmpi(BasisFunc,'Harmonic ODE')
        lsqfunc= @(Coef)fModel(Coef)-yData;
    else
        lsqfunc= @(Coef)fModel(Coef,xData)-yData;
    end
    init=rand(1,numCoef);
    lsqBeta=lsqnonlin(lsqfunc,init);
    if strcmpi(Solver,'Lsqnonlin')
        Beta=lsqBeta;
    end
end
%% fminsearch
if strcmpi(Solver,'fminsearch') || strcmpi(Solver,'All')
    if strcmpi(BasisFunc,'Harmonic ODE')
        fR= @(Coef)fModel(Coef)-yData;
    else
        fR=@(Coef)fModel(Coef,xData)-yData;
    end
    fminfunc= @(Coef)fR(Coef)'*fR(Coef);
    options = optimset('TolFun',1e-4,'TolX',1e-4);
    init=rand(1,numCoef);
    fminBeta=fminsearch(fminfunc,init,options);
    if strcmpi(Solver,'fminsearch')
        Beta=fminBeta;
    end
end
%% Display Results
switch Solver
    case 'all'
        DispLinBeta=sprintf('%.2f,',linBeta);
        DisplsqBeta=sprintf('%.2f,',lsqBeta);
        DispfminBeta=sprintf('%.2f,',fminBeta);
        DispMessage=sprintf('Model Coefficients for each Solver\nLinear: %s\nlsqnonlin: %s\nfminsearch: %s',DispLinBeta,DisplsqBeta,DispfminBeta);
        disp(DispMessage)
    case 'fminsearch'
        DispfminBeta=sprintf('%f,',fminBeta);
        DispMessage=sprintf('Model Coefficients for fminsearch\n%s',DispfminBeta);
        disp(DispMessage)
    case 'lsqnonlin'
        DisplsqBeta=sprintf('%f,',lsqBeta);
        DispMessage=sprintf('Model Coefficients for lsqnonlin\n%s',DisplsqBeta);
        disp(DispMessage)
    otherwise
        fprintf('Error!! Solver not recognized')
        keyboard
end


if strcmpi(PlotBool,'yes') || strcmpi(PlotBool,'no')
    if strcmpi(Solver,'All')
        %Plotting Commands for if plotting data from all solvers
    else
        
    end
    if strcmpi(BasisFunc,'Harmonic ODE')
        [xPlot,yPlot]=HarmonicODE(Beta,tspan);
    else
        xPlot=linspace(min(xData)-(min(xData)+max(xData))/10,...
        max(xData)+(min(xData)+max(xData))/10);
        yPlot=fModel(Beta,xPlot);
    end
    plot(xData,yData,'r*','MarkerSize',10)
    hold on 
    plot(xPlot,yPlot)
    switch BasisFunc
        case '0 Order Poly'
            LegendFunc=sprintf('y=%.2f',Beta);
        case '1st Order Poly'
            LegendFunc=sprintf('y=%.2f%+.2fx',Beta);
        case '2nd Order Poly'
            LegendFunc=sprintf('y=%.2f%+.2fx%+.2fx^2',Beta);
        case '3rd Order Poly'
            LegendFunc=sprintf('y=%.2f%+.2fx%+.2fx^2%+.2fx^3',Beta);
        case 'Mixed Exponential'
            LegendFunc=sprintf('y=%.2f%+.2fe^{%.2fx%+.2fx^2}',Beta);
        case 'Harmonic ODE'
            LegendFunc=sprintf('y''''%+.2fy''%+.2fy=0',Beta(3),Beta(4));
        case 'Self Enter'
        otherwise
            fprintf('Error!! BasisFunc not recognized')
            keyboard
    end
    legend('Data',LegendFunc)
    xlabel('x')
    ylabel('y')
    PlotTitle=sprintf('%s Function Fit for %s',BasisFunc,DataType);
    title(PlotTitle)
    
end
%%
function [x,y]=HarmonicODE(Coef,tspan,varargin)
    y0=[Coef(1),Coef(2)];
    dydt = @(t,y)[y(2); -Coef(3).*y(2) - Coef(4).*y(1)];
    [x,y]=ode45(dydt,tspan,y0);
    if size(varargin)~=0
       NumPoints=varargin{1};
       index=floor((1:NumPoints)*length(x)/20); 
       x=x(index);
       y=y(index,1);
    end
end
    



