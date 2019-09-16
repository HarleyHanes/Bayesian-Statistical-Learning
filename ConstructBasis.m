function [Basis] = ConstructBasis(Basis)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
switch Basis.Func
    case '0 Order Poly'
        fModel = @(Coef)Coef(1);
        numCoef=1;
    case '1st Order Poly'
        fModel = @(Coef)Coef(1) + Coef(2)*x;
        numCoef=2;
    case '2nd Order Poly'
        fModel = @(Coef)Coef(1) + Coef(2)*x +Coef(3)*x.^2;
        numCoef=3;
    case '3rd Order Poly'
        fModel = @(Coef,x)Coef(1) + Coef(2)*x +Coef(3)*x.^2 +Coef(4)*x.^3;
        numCoef=4;
    case 'Mixed Exponential'
        fModel = @(Coef,x)Coef(1) + Coef(2) * exp(Coef(3).*x + Coef(4)*x.^2);
        numCoef=4;
    case 'First Order ODE'
        fModel=@(Coef,tspan)FirstOrderODE(Coef,tspan);
        numCoef=2;
    case 'Harmonic ODE'
        fModel=@(Coef,tspan)HarmonicODE(Coef,tspan);
        numCoef=4;
    case 'Self Enter'
        fModel=Basis.fDefine;
    otherwise
        fprintf('Error!! BasisFunc not recognized')
        keyboard
end
Basis.fModel=fModel;
Basis.numCoef=numCoef;
end

