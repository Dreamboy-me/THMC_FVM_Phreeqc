function [muW1,rhoW1] = fluidProperties(G,pNew,TNew,dynamic)
if dynamic == 0
    muW1 = 1*centi*poise.*ones(G.cells.num,1);
    rhoW1 = 1000*kilogram/meter^3.*ones(G.cells.num,1);
else
    % Inspect fluid model
    % the fluid density and viscosity as a function of pressure and
    % temperature. 
    setAxProps = @(ax) set(ax, 'View'              , [-140,20]    , ...
        'PlotBoxAspectRatio', [2,2,1]      , ...
        'Projection'        , 'Perspective', ...
        'Box'               , 'on'         , ...
        'XLimSpec'          , 'tight'      , ...
        'YLimSpec'          , 'tight'      , ...
        'ZLimSpec'          , 'tight'      );

    fluid = initSimpleADIFluid('mu', 0.001, 'rho', 1000, 'phases', 'W');
    fluid = addThermalFluidProps(fluid, 'Cp'     , 4.0e3, ...
        'lambdaF', 0.5  , ...
        'useEOS' , true );

    K0   = 273.15*Kelvin;
    pMin = 1e6*Pascal;      % Minimum pressure
    pMax = 200e6*Pascal;    % Maximum pressure
    TMin = K0;              % Minimum temperature
    TMax = K0 + 275*Kelvin; % Maximum temperature
    n    = G.cells.num;

    % Get pressure/temperature grid
    p1    = pNew;
    T1    = K0+TNew;
    rhoW1 = fluid.rhoW(p1, T1);
    muW1 = fluid.muW(p1, T1);
end
end
