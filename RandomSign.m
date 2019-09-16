function [Correlated,z] = RandomSign(R,Signif, displayResults)
%RandomSign Residual Analysis test for Randomness of signs
%   Authors: Harley Hanes, Elliot Hill, Helen Weierbach
%   Inputs: A vector "R" of the residuals 
%           A scalar "Signif" identifying the minimum z score for the
%            signs to be denoted as correlated
%           A boolean "displayResults" denoting whether a message detailing
%            z score and and pass/ fail should be printed
%   Outputs:A boolean "Correlated" giving whether the test passed or failed 
%            at the given significance level
%           A scalar "z" of the statisitical z-score of likeliness on a 
%            standard normal distribution
%   References: Hansen C., Pereyra V., Scherer G. 2012. Least Squares Data
%               Fitting with Application. Johns Hopkins University Press


%Count culmulative signs, and changes
m=length(R);
nPositive=sum(R>0);
nNegative=sum(R<0);
u = sum(abs(diff(sign(R))))/2+1;          % Find the number of sign shifts 
                                          % and add 1 to get the number of 
                                          % runs

%compute mean
mean=(2*nPositive*nNegative)/m+1;          %See 1.4.1 from p.14 of Hensen
%compute standard deviation
stDeviation=sqrt((mean-1)*(mean-2)/(m-1)); %See 1.4.1 from p.14 of Hensen
%compute score
z=abs(u-mean)/stDeviation;

%Check if significant
if z < Signif
        Correlated=0;
    else
        Correlated=1;
end

%Display Results
if displayResults == true
    disp(['z = ', num2str(z)]);
    if Correlated==0
        fprintf(['z is less than significance level of %.2f so the '...
            'signs are not \nsignificantly correlated.\n'],Signif);
    else
        fprintf(['z is greater than significance level of %.2f so the '...
            'signs are \nsignificantly correlated.\n'],Signif);
    end
end

end

