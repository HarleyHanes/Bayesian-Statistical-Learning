%InvestigateSVD
%Harley Hanes
%Test code for investigating properties of SVD decomp
clc;clear;

%% Generating A
mRow=3;
nCol=5;
rowsingular=1;
colsingular=1;
A=rand(mRow,nCol);
if rowsingular==1
    A(end,:)=A(end-1,:);
end
if colsingular==1
    A(:,end)=A(:,end-1);
end
disp('A='), disp(A)

%% Assessing SVD
disp('SVD Properties')
[U,D,V]=svd(A);  %A=UDV' not UDV
ASVD=U*D*V';
SVDerr=norm(A-ASVD);
disp('SVD Error:'),disp(SVDerr)

%Assessing rank
    Arank=rank(A); disp('Rank of A:'),disp(Arank)
    Urank=rank(U); disp('Rank of U:'),disp(Urank)
    Drank=rank(D); disp('Rank of D:'),disp(Drank)
    Vrank=rank(V); disp('Rank of V:'),disp(Vrank)
    %--rank(U)=mRow, rank(V)=nCol, rank(D)=rank(A)

    %How does At compare to A?
    At=A';disp('A transpose:'),disp(At)
    [Ut,Dt,Vt]=svd(At);
    Atrank=rank(At); disp('Rank of At:'),disp(Atrank)
    Utrank=rank(Ut); disp('Rank of Ut:'),disp(Utrank)
    Dtrank=rank(Dt); disp('Rank of Dt:'),disp(Dtrank)
    Vtrank=rank(Vt); disp('Rank of Vt:'),disp(Vtrank)
    %--A, A', D, and Dt have same rank as expected but the ranks of U and V
    %   switch for Ut and Vt which makes sense since the column space will switch
    %   will the row space.
%% Assessing MPP
disp('MPP Properties')
MPP=A'*A;disp('MPP:'),disp(MPP)
%How does A'*A compare to A*A'?
AAt=A*A';disp('AAt:'),disp(AAt)

    %Rank
    MPPrank=rank(MPP);disp('Rank of MPP:'),disp(MPPrank)
    AAtrank=rank(AAt);disp('Rank of AAt:'),disp(AAtrank)
    %--They have the same rank as expected
    
    %Column and row sums
    ColSumMPP=sum(MPP);disp('Column Sum of MPP:'),disp(ColSumMPP)
    RowSumMPP=sum(MPP');disp('Row Sum of MPP:'),disp(RowSumMPP)
    ColSumAAt=sum(AAt);disp('Column Sum of AAt:'),disp(ColSumAAt)
    RowSumAAt=sum(AAt');disp('Row Sum of AAt:'),disp(RowSumAAt)
    %--They have the same column and row sums
%% How do eigenvalues and vectors of MPP and AAt compare to SVD(A)?
    [VMPP,DMPP,WMPP]=eig(MPP); %--Symmetric so VMPP=WMPP
    [VAAt,DAAt,WAAt]=eig(AAt); %--Symmetric so VAAt=WAAt
    %compare matrices eigenvectors matrices of equal size
    %--For and SVD matrix Z and eigenvector matrix of equal size V from MPP
    %  or AAt, Z is almost left inverse of V' but only if Z is full rank
    %  and even then the signs are off and its on the off diagonal
    %--If Z isn't full rank then the property holds but only for the first
    %  r rows where r=rank(Z)
        order={'VMPP','VAAt'};
        Vlist={VMPP,VAAt};
        for i=1:2
            Veig=Vlist{i};
            if size(Veig)==size(U)
                fprintf('Rank of %s:',order{i}),disp(rank(Veig))
                disp('Rank of U:'),disp(rank(U))
                fprintf('%s-U:',order{i}),disp(Veig-U)
                fprintf('%s*U'':',order{i}),disp(Veig*U')
                fprintf('%s''*U:',order{i}),disp(VAAt'*U)
            end
            if size(Veig)==size(V)
                fprintf('Rank of %s:',order{i}),disp(rank(Veig))
                disp('Rank of V:'),disp(rank(V))
                fprintf('%s-V:',order{i}),disp(Veig-V)
                fprintf('%s*V'':',order{i}),disp(Veig*V')
                fprintf('%s''*V:',order{i}),disp(V'*Veig)
            end
        end









