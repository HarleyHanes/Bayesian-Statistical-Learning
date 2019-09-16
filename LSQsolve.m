function [Results] = LSQsolve(Data,Basis,Solver)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
fModel=Basis.fModel;
numCoef=Basis.numCoef;
%% Linear Solver
if strcmpi(Solver,'Linear') || strcmpi(Solver,'All')
    rho=@(x)fModel(ones(1,numCoef),x);
    X=rho(Data.x);
    linBeta=((X'*X)\(X'*Data.y))';
    if size(linBeta)== size(rho(Data.x(1)))
        Results.linBeta=linBeta';
    end
    if strcmpi(Solver,'linear')
        Results.Beta=linBeta;
    end
end

%% Lsqnonlin
if strcmpi(Solver,'Lsqnonlin') || strcmpi(Solver,'All')
    lsqfunc= @(Coef)fModel(Coef,Data.x)-Data.y;
    init=rand(1,numCoef);
    Results.lsqBeta=lsqnonlin(lsqfunc,init);
    if strcmpi(Solver,'Lsqnonlin')
        Results.Beta=Results.lsqBeta;
    end
end
%% fminsearch
if strcmpi(Solver,'fminsearch') || strcmpi(Solver,'All')
    fR=@(Coef)fModel(Coef,Data.x)-Data.y;
    fminfunc= @(Coef)fR(Coef)'*fR(Coef);
    options = optimset('TolFun',1e-4,'TolX',1e-4);
    init=rand(1,numCoef);
    Results.fminBeta=fminsearch(fminfunc,init,options);
    if strcmpi(Solver,'fminsearch')
        Results.Beta=Results.fminBeta;
    end
end
end

