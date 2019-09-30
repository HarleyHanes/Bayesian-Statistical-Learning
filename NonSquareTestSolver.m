%%Solving non-square systems
%-Harley Hanes, Fall 2019
%Testing SVD, Morse-Penrose Psuedoinverse (MPP), and Ridge/ Tickinoff
%methods for solving nonsquare systems. 
%% Results Notes
%Worried There may be an issue with my functions because SVD and MPP are
%failing for full rank matrices
% m>n
% 1)A\b, SVD and MPP methods working for full rank matrix. Ridge Fails
% 2)A\b has err of 3.92 (in A*x), SVD has error 33.63, MPP has error of 4.58.
%                          Ridge has err of 8.45
%                          MPP gives very large value (O(10^14)) in x
% m<n
% 1)Only A\b works for full rank. SVD has error in A*x of 15.9, MPP of 108
%                                 Ridge of 7.999
%                                 SVD has erro in x of 6.0753, MPP of 43.7
%                                 Ridge has 4.48
% 2)A\b has err in A*x of .7071. SVD has error on order of 10^14, MPP has
%                                error of 137, Ridge has 7.61
% m=n
% 1) Only A\b works for full rank. SVD and MPP have error of 1.38 (same
%                                   results), Ridge has 1.73
% 2) Only A\b works for rank deff. SVD has error of 10.9, MPP of 1.58
%                                  Ridge has err of 1.73 (gives all ones)
                                  

%%
clear;clc;
%Problem Set Up
%m>n
% A=[1 2; 3 4; 5 6]; %Declare Full rank matrix
% A=[1 2; 3 4; 5 6]; A(:,1)=A(:,2);  %Declare rank defecient matrix
% b=[7;8;9];                           %Declare Output
%m<n
% A=[1 2 3; 4 5 6];
% A=[1 2 3; 1 2 3];
% b=[7;8];
% m=n
% A=[1 2 0; 0 3 4; 5 0 6];
A=[1 2 3; 4 5 6; 7 8 9];
b=[1; 1; 1];
%SVD
ASVD=SVDinv(A);
xSVD=ASVD*b;
%MPP
MPP=A'*A;
AMPP=SVDinv(MPP);
xMPP=AMPP*A'*b;
%Ridge/Tickinoff
xRidge=ridge(b,A,1);
%Results Comparison
x=[A\b xSVD xMPP xRidge];
xerr=b-A*x;
xerr(end+1,:)=sqrt(sum(xerr.^2));
%Functions
%SVD
function Ainv=SVDinv(A)
    [U,D,Vt]=svd(A);
    V=Vt';
    Dinv=zeros(size(D'));
    for i= 1:min(size(D))
        if D(i,i)>0
            Dinv(i,i)= 1/ D(i,i);
        end
    end
    Ainv=V*Dinv*U';
end