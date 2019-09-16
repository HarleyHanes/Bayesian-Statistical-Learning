function Position=FirstOrderODE(Coef,tspan,varargin)
    y0=Coef(1);
    dydt=@(t,y)Coef(2)*y;
    [x,Position]=ode45(dydt,tspan,y0);
     x;
end

