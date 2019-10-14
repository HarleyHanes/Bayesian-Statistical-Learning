function x = ScaleSolve(A,b,Abnorm)
%ScaleSolve TestCode for Solving a linear system with scaling
%   Inputs: A equation matrix, b condition vector, norm to scale A and b
%   Outputs: x solution vector, unscaled
%   Solving Method: linsolve
%   Authors: Harley Hanes

%% scale then ridge
 bScaled=b/norm(b,Abnorm);         %Scale b to unit
 Atilda=A/norm(b,Abnorm);          %Account for b scaling in A
 [~,phisize]=size(A);
 phi=NaN(1,phisize)';             %Set up A Scaling matrix
 for i=1:length(phi)
     for j=1:length(phi)
         phi(j)=1/norm(Atilda(:,j),Abnorm); 
     end
 end
     AScaled=Atilda.*phi';
     
     
     
     %% Check A and b scalings
     for j=1:length(A)
         if norm(AScaled(:,j),Abnorm) < 1-10^(-6) || norm(AScaled(:,j),Abnorm) > 1+10^(-6)
             warning('Atilda Scaling Unsuccesful')
             keyboard
         end
     end
     if norm(bScaled,Abnorm) < 1-10^(-6) || norm(bScaled,Abnorm) > 1+10^(-6)
        warning('b Scaling Unsuccessful')
        keyboard
     end
     
     %Scale x
     xScaled=linsolve(AScaled,bScaled);
     Axnorm=norm(A*xScaled);
     %xScaled=xScaled*norm(b,Abnorm);
     
     %Check x Scalings
     if Axnorm < 1-10^(-6) || Axnorm > 1+10^(-6)
        warning('x Scaling Unsuccssful')
        keyboard
     end
     
     x=xScaled.*phi;
end