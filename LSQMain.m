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
% 1)All Solvers, DataTypes, and BasisFunctions working and plotting
% correctly
%===================================================================
clear; close all;clc;
%% Master Control
Data.Type='Harmonic ODE'; %Options: 'Class Data', 'Normal Scatter
                          %         'Straight Line', 'Runge
                          %         Function','First Order ODE','Harmonic
                          %         ODE','Load Data'
   
   if strcmpi(Data.Type,'Load Data')
       Data.filename='data.csv';%Include full path name if not in path
       Data.DataOrient='column';    %whether data vectors are column or row vectors
       Data.xPoint=1;           %which row/ column xData is in
       Data.yPoint=2;           %which row/ column yData is in
   else
       Data.xMin=0;
       Data.xMax=10;
       Data.numPoints=20;
       Data.Coef=[0 4 1 2];
   end
Basis.Func='3rd Order Poly'; %Opions '0, 1st, 2nd, 3rd Order Poly'
                            %        'Mixed Exponential','First Order ODE'
                            %        'Harmonic ODE','Define Own'
    if strcmpi(Basis.Func,'Define Own')
        Basis.fDefine=@(Coef,x) Coef(1).*x;
        Basis.numCoef=1;
    end
Solver='fminsearch'; %Options: Linear, lsqnonlin, fminsearch, all
PlotBool='yes';


%% Data and basis function entry
Data=ConstructData(Data);
Basis=ConstructBasis(Basis);

%Note: Data vectors need to be column vectors for fminsearch and linear

%% Solver
Results=LSQsolve(Data,Basis,Solver);


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

    xPlot=linspace(min(xData),max(xData));
    yPlot=fModel(Beta,xPlot);
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
        case 'First Order ODE'
            LegendFunc=sprintf('y_0=%.2f, y''=%.2fy',Beta);
        case 'Harmonic ODE'
            LegendFunc=sprintf('y_0=%.2f,y_0''=%.2f\ny''''%+.2fy''%+.2fy=0',Beta);
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
%% ODE functions  



