%%DifferentialAlgebraTestCode
%-Harley Hanes, Fall 2019
% Code testing

clear;clc;close all
A=[-4 1; 1 2];
B=[1 0];
n=2:120;
for i=1:length(n)
    t=linspace(0,1,n(i));
    y=linspace(0,2,n(i))';
    x0(:,i)=DiffAlgebra(A,B,t,y);
end

plot(n,x0)

function x0=DiffAlgebra(A,B,t,y)
    if length(A)~=length(B)
        fprintf('Error!! A and B need the same number of columns.')
        keyboard
    end
    [U,D]=eig(A);
    S=NaN(length(t),length(A));
    for i=1:length(t)
        S(i,:)=B*(U*exp(D*t(i))/U);
    end
    x0=S\y;

end

