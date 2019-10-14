function x = RidgeScale(A,b,Abnorm,lambda)
%RidgeScale TestCode for Solving a Ridge decomposition with scaling
%   Inputs: A equation matrix, b condition vector, norm to scale A and b,
%   lambda for use in ridge
%   Outputs: x solution vector, unscaled
%   Authors: Harley Hanes

%% scale then ridge
% bScaled=b/norm(b,Abnorm);         %Scale b to unit
% Atilda=A/norm(b,Abnorm);          %Account for b scaling in A
% phi=NaN(length(A),length(A'));             %Set up A Scaling matrix
% for i=1:length(phi')
%     for j=1:length(phi)
%         phi(:,j)=1/abs(Atilda(j,i)); 
%     end
% end
% AScaled=Atilda*phi;
%     AL=AScaled+eye(length(AScaled))*lambda;   %Apply Ridge (no longer Unit 
%                                               %though)
%     ARidge=AL'*AL;
%     bRidge=AL'*bScaled;
%     ARidgeInv=SVDinv(ARidge);
%     xScaled=ARidgeInv*bRidge;
% x=phi*xScaled;

%% Ridge then Scale
%Make Ridge matrices
AL=A+eye(length(AScaled))*lambda;
ARidge=AL'*AL;
bRidge=Al'*b;
%Scale
bScaled=bRidge/norm(b,Abnorm);         %Scale b to unit
Atilda=ARidge/norm(b,Abnorm);          %Account for b scaling in A
phi=NaN(length(Atilda));             %Set up A Scaling matrix
for i=1:length(phi)
    for j=1:length(phi)
        phi(:,j)=1/abs(Atilda(j,i)); 
    end
end
AScaled=Atilda*phi;
xScaled=SVDinv(AScaled)*bScaled;
x=phi*xscaled;

end

