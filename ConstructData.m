function Data = ConstructData(Data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
DataType=Data.Type;
if strcmpi(DataType,'Load Data')
    filename=Data.filename;
    DataOrient=Data.DataOrient;
    xPoint=Data.xPoint;          
    yPoint=Data.yPoint;
else
    xMin=Data.xMin;
    xMax=Data.xMax;
    numPoints=Data.numPoints;
    xData=linspace(xMin,xMax,numPoints)';
    Coef=Data.Coef;
end
%filename and DataOrient declared below so they are not required for other
%Data type entries
switch DataType
    case 'Normal Scatter'
        yData=Coef(1)+randn(length(xData),1)*Coef(2);
    case 'Straight Line'
        yData=Coef(1)*ones(length(xData),1)+Coef(2);
    case 'Runge Function'
        a=Coef(1);
        b=Coef(2);
        c=Coef(3);
        TrueFunc=@(x) a./(b+c.*x.^2);
        yData=TrueFunc(xData);
    case 'First Order ODE'
        yData=FirstOrderODE(Coef,xData);
    case 'Harmonic ODE'
        yData=HarmonicODE(Coef,xData);
    case 'Load Data'
        DataTemp=load(filename);
        switch DataOrient
            case 'column'
                xData=DataTemp(:,xPoint);
                yData=DataTemp(:,yPoint);
            case 'row'
                xData=DataTemp(xPoint,:);
                yData=DataTemp(yPoint,:);
            otherwise
                fprintf(['Error!! Data orientation not recognized. Enter as'...
                    'row or column'])
                keyboard
        end
    otherwise 
        fprintf('Error!! DataType not recognized')
        keyboard
end
Data.x=xData;
Data.y=yData;
end

