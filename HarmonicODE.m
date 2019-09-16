function Position=HarmonicODE(Coef,tspan)
    y0=[Coef(1),Coef(2)];
    dydt = @(t,y)[y(2); -Coef(3).*y(2) - Coef(4).*y(1)];
    [x,y]=ode45(dydt,tspan,y0);
    Position=y(:,1);
    x;
end
