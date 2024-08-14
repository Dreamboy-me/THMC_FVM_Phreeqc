%% Single-Phase Problem Demonstrating the sequential algorithm for THMC coupled modeling Model
% In this example, we will demonstrate how to model the THMC coupled 
% process for EGS using a unified finite volume method (FVM) computational 
% framework, which employs a sequential and one-way coupling mechanism to 
% model the hydraulic, thermal, chemical, and mechanical processes. 
% Chemical reactions within the solute transport module are 
% simulated using PHREEQC through a dedicated interface.
% phreeqc interface is ok 
close all; clear;
pause(2)
clc;
pause(3)

startsub;
pause(2)
% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add incomp;          % Incompressible fluid models
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture lines
% Construct a Cartesian grid comprising 151-by-151 cells, where each cell has
% dimension 0.6623-by-0.6623 m^2. 

celldim = [151, 151]; % number of cells--> not nodes
physdim = [100, 100];

nx = celldim(1); ny = celldim(2);
L = physdim(1);  D = physdim(2); 
dx = L/nx;   dy = D/ny;

GG = cartGrid(celldim, physdim);
Gmfr = computeGeometry(GG);

%% fracture coords
%  fractures lines in [x1 y1 x2 y2] format.
fl  = [1.5, 1.2, 98.5, 98.2;
      ...
      ]; % hydaulic fracture and natural fracture  

%% Process fracture lines
% Using the input fracture lines, identify independent fracture networks
% comprising of connected lines. 
% The matrix grid and fracture lines are plotted.

dispif(mrstVerbose, 'Processing user input...\n\n');
[Gmfr,fracture] = processFracture2D(Gmfr,fl);
fracture.aperture = 1/1000; % Fracture aperture

figure;
plotFractureLines(Gmfr,fracture);
axis equal tight;
box on
dispif(true, 'rock grids and fracture lines generation finished\n\n');

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model).

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
Gmfr = CIcalculator2D(Gmfr,fracture);
[Gmfr,F,fracture] = gridFracture2D(Gmfr,fracture,'assemblyType',3); % only this method

clf; plotFractureNodes2D(Gmfr,F,fracture,'linewidth',0.5,'markersize',2.0,'shownumbering',0);
axis equal tight; box on
dispif(true, 'CI finished\n\n'); 

%% Set rock properties in fracture and matrix
% Set the permeability (K)  in the matrix and in the fractures.

GC = Gmfr; GT = Gmfr;
dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');

% non_uniform of permeability and porosity
p = gaussianField(Gmfr.cartDims, [0.035 0.045], [11 3]); 
K = p.^3.*(1e-5)^2./(5*2.315^2*(1-p).^2);
Gmfr.rock.poro = p(:);   
Gmfr.rock.perm = K(:);   

figure
plotCellData(Gmfr, Gmfr.rock.perm,'EdgeColor','none');
colormap jet(25)
view(0, 90); colorbar; axis equal tight
figure
plotCellData(Gmfr, Gmfr.rock.poro,'EdgeColor','none');
colormap jet(25)
view(0, 90); colorbar; axis equal tight

poro_frac = 0.8;
K_frac = 1.0e-11/darcy;  
Gmfr = makeRockFrac(Gmfr, K_frac, 'permtype','homogeneous','porosity', poro_frac); 

cond_fluid = 0.5; % fluid conductivity coefficent
cond_solid = 3.0; % solid conductivity coefficent
GT.rock.poro = Gmfr.rock.poro; %  Matrix_porosity 

GT.rock.perm = GT.rock.poro*cond_fluid+(1-GT.rock.poro)* cond_solid;  % Matrix efficient_λ
KT_frac = poro_frac*cond_fluid/darcy +(1-poro_frac)*cond_solid/darcy; % Fracture efficient_λ
GT = makeRockFrac(GT, KT_frac, 'permtype','homogeneous','porosity', poro_frac);

GC.rock.poro = Gmfr.rock.poro; % Matrix_porosity 
Diff = 1e-9; 
GC.rock.perm = GC.rock.poro * Diff; %  φ_Matrix * D
KC_frac = poro_frac * Diff/darcy;   %  φ_Fracture * D
GC = makeRockFrac(GC, KC_frac, 'permtype','homogeneous','porosity', poro_frac); 
dispif(true, 'Porosity and permeability finished\n\n');

%% Define fracture connections as NNC and compute the transmissibilities
% In this section, use the function defineNNCandTrans to combine the
% fracture and matrix grid structures into a single grid structure. In
% addition to that, assign a 'non-neighbouring connection (NNC)' status
% to every fracture-matrix connection. 
% To compute the flux between these elements, compute a
% transmissibility for each NNC using the CI's computed earlier. 
% Vector T contains the transmissibility for each face 
% in the combined grid and each NNC.

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
[Gmfr,T] = defineNNCandTrans(Gmfr,F,fracture);
[GT,TT]  = defineNNCandTrans(GT,F,fracture);  % added-->k insteaded by λ
[GC,TC]  = defineNNCandTrans(GC,F,fracture);  % added-->k insteaded by D

%% Add BC
% Set boundary condition 
%
% Pressure
bc = [];

% wellbore
inj = 1;
prod = Gmfr.cartDims(1) * Gmfr.cartDims(2);
wellRadius = 0.1;
W = addWell([], Gmfr.Matrix, Gmfr.Matrix.rock, inj,'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', 220*meter^3/day, 'Radius', wellRadius); 
W = addWell(W, Gmfr.Matrix, Gmfr.Matrix.rock, prod, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 47.7*barsa, 'Radius', wellRadius, 'Type_well', 'production');

% Temperature
bcT = [];

% Concentration
bcC = [];

%% Define fluid properties
% Define a single fluid of viscosity and density.

fluid  = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3); 
fluidT = initSingleFluid('mu', 1.0, 'rho', 2623*kilogram/meter^3); % added-->T,unsteady,mu insteaded by fluid ρ
fluidC = initSingleFluid('mu', 1.0, 'rho', 2623*kilogram/meter^3); % added-->C,unsteady
state_P  = initResSol(Gmfr, 30.0*barsa);
state_T  = initResSol(GT, 236.0); % -->T 
state_C  = initResSol(GC, 0.0);  % -->C

thermophysical.cpf = 4000;
thermophysical.cpr = 980;
Cpm = thermophysical.cpr ; % Matrix,Fracture,thermalphysical_cp
Cpf = thermophysical.cpf;  % Fluid_thermalphysical_cp
thermophysical.rhor = 2623;
thermophysical.conductivityf = cond_fluid;
thermophysical.conductivityr = cond_solid;
betaT = 5.0e-6;
Ct = 4.4e-5/barsa; 
Tinj = 40; 
Cinj = 1.0e-7; 

%% Time step added

numSteps = 40;               % number of time-steps
totTime  = numSteps*360*day; % total simulation time 
tol      = 1e-5;                % Newton tolerance
maxits   = 100;                 % max number of Newton its
rampSteps = numSteps; ddt = totTime / numSteps; % dt_list = repmat(dt,numSteps,1); 
dt_list = rampupTimesteps(totTime, ddt, rampSteps); % simple rampup steps
totNumdt = numel(dt_list); 
time_list = cumsum(convertTo(dt_list,year));

%% Loop initial value added 
% Initial conditions are adjusted in response to real-time changes in the simulation conditions

prNew    = 7.0e6.*ones(Gmfr.Matrix.cells.num,1); 
pfNew    = 7.0e6.*ones(Gmfr.cells.num-Gmfr.Matrix.cells.num,1);
pinj     = 7.0e6;

TrNew    = 236.*ones(Gmfr.Matrix.cells.num,1);
TfNew    = 236.*ones(Gmfr.cells.num-Gmfr.Matrix.cells.num,1);

CrNew_SiO2  = 0.0.*ones(GC.Matrix.cells.num,1);
CfNew_SiO2  = 0.0.*ones(GC.cells.num-GC.Matrix.cells.num,1);

CrNew_K = 0.0e-5.*ones(GC.Matrix.cells.num,1);
CfNew_K = 0.0e-5.*ones(GC.cells.num-GC.Matrix.cells.num,1);

CrNew_Na = 0.0.*ones(GC.Matrix.cells.num,1);
CfNew_Na = 0.0.*ones(GC.cells.num-GC.Matrix.cells.num,1);

CrNew_Al = 0.0e-5.*ones(GC.Matrix.cells.num,1);
CfNew_Al = 0.0e-5.*ones(GC.cells.num-GC.Matrix.cells.num,1);

CrNew_Mg = 0.0.*ones(GC.Matrix.cells.num,1);
CfNew_Mg = 0.0.*ones(GC.cells.num-GC.Matrix.cells.num,1);

CrNew_H = 0.0e-7.*ones(GC.Matrix.cells.num,1);
CfNew_H = 0.0e-7.*ones(GC.cells.num-GC.Matrix.cells.num,1);

pHr = 7.0*ones(GC.Matrix.cells.num,1);
pHf = 7.0*ones(GC.cells.num-GC.Matrix.cells.num,1);

poroC = Gmfr.rock.poro(:);
permC = Gmfr.rock.perm(:);

Tri = TrNew;
Tfi = TfNew;
% The volume fraction of each component within the mineral
quartz = 0.341; kspar = 0.15; albite = 0.321; phlogopite = 0.168; clinochlore = 0.02; 
GC = chemical_initial(GC, quartz, kspar, albite, phlogopite, clinochlore);     % 矿物的体积、石英等密度和离子浓度

quartzC = GC.chemical.quartz(:); 
ksparC = GC.chemical.kspar(:); 
albiteC = GC.chemical.albite(:);  
phlogopiteC = GC.chemical.phlogopite(:); 

%% Compute process

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
sol = repmat(struct('time',[], 'prNew',[],'pfNew',[],'pinj',[],'Fluxinj',[],'Fluxout',[],'TrNew',[],'TfNew',[],'Pwinj',[],'Pwout',[],'Tout',[],'CrNew_SiO2',[],'CfNew_SiO2',[],'CrNew_K',[],'CfNew_K',[],'CrNew_Na',[],'CfNew_Na',[],'CrNew_Al',[],'CfNew_Al',[],'CrNew_Mg',[],'CfNew_Mg',[],'CrNew_H',[],'CfNew_H',[],'pHr',[],'pHf',[],'quartzC',[],'ksparC',[],'albiteC',[],'phlogopiteC',[],'poroC',[],'permC',[],'u_x',[],'u_y',[]),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'prNew', prNew,'pfNew', pfNew,'pinj', pinj,'Fluxinj',0,'Fluxout',0,'TrNew',TrNew,'TfNew',TfNew,'Pwinj',0,'Pwout',0,'Tout',0,'CrNew_SiO2',CrNew_SiO2,'CfNew_SiO2',CfNew_SiO2,'CrNew_K',CrNew_K,'CfNew_K',CfNew_K,'CrNew_Na',CrNew_Na,'CfNew_Na',CfNew_Na,'CrNew_Al',CrNew_Al,'CfNew_Al',CfNew_Al,'CrNew_Mg',CrNew_Mg,'CfNew_Mg',CfNew_Mg,'CrNew_H',CrNew_H,'CfNew_H',CfNew_H,'pHr',pHr,'pHf',pHf,'quartzC',quartzC,'ksparC',ksparC,'albiteC',albiteC,'phlogopiteC',phlogopiteC,'poroC',poroC,'permC',permC,'u_x',0,'u_y',0);

states = deal(cell(totNumdt,1));

t = 0; step = 0; 
hwb = waitbar(0,'EGS THMC Simulation...miracle will happen');
while (t<totTime) && (step < totNumdt)
tic;
% while t<totTime
    step = step + 1;
    dt = dt_list(step); % for unconstan time step 
    dtc    = dt_list(step)/iterc; % for unconstant time step

    t = t + dt;
    fprintf('\nTime step %d: Time %.2f -> %.2f years\n',step, convertTo(t - dt, year), convertTo(t, year)); 
    % Newton loop
    resNorm =  1e99;    
    prOld = prNew;
    pfOld = pfNew;
    
    TrOld = TrNew;
    TfOld = TfNew;    
   
    pNew = [prNew;pfNew];
    TNew = [TrNew;TfNew];
    
    [viscosity,density] = fluidProperties(Gmfr,pNew,TNew,1); % Physical properties vary with changes in temperature and pressure.
    [prfv,srfv]   = poreVolumeDynamic(Gmfr); % The variation of pore volume with pressure
    nit = 0;
    while (resNorm > tol) && (nit < maxits)
        Pr = pNew;
        Tr = TNew;
        
        % Fluid flow
        stateP = incompTPFAFVM(state_P, Gmfr, T, fluid,viscosity,density,prfv,Ct,prOld,pfOld,dt,'bc',bc,'Wells', W, 'gravity',gravity,'use_trans',true);
        prNew = stateP.prNew;
        pfNew = stateP.pfNew;
        pNew = [prNew;pfNew]; 
        pinj = stateP.wellSol(1).pressure;
        Fluxinj = stateP.wellSol(1).flux;
        Fluxout = stateP.wellSol(2).flux;
        velocity  = stateP.velocity;
        velocityC = stateP.velocityC;
        
        % Heat transfer
        stateT = incompTPFATWell(state_T, GT, TT, fluidT, fluid, velocity, prNew, TrOld, TfOld, density, viscosity, thermophysical, Cpm, Cpf, prfv, srfv, nx, ny, dt, Tinj, 'bc', bcT, 'Wells', W, 'use_trans', true);
        TrNew = stateT.TrNew;
        TfNew = stateT.TfNew; 
        TNew  = [TrNew;TfNew]; 

        Pwinj = stateT.wellSol(1).thermalpower;
        Tout  = stateT.wellSol(2).temperature;
        Pwout = stateT.wellSol(2).thermalpower;

        % Reactive solute transport  
        for i=1:iterc            
            GC = incompTPFACWell(state_C, GC, F, fracture, fluidC, velocityC, prNew,  pfNew, TrNew, TfNew, density, viscosity, Diff, dtc, Cinj, bcC, W);
        end  
        %
        CrNew_SiO2 = GC.chemical.ion.SiO2(1:GC.Matrix.cells.num);
        CfNew_SiO2 = GC.chemical.ion.SiO2(GC.Matrix.cells.num+1:end);
        
        CrNew_K = GC.chemical.ion.K(1:GC.Matrix.cells.num);
        CfNew_K = GC.chemical.ion.K(GC.Matrix.cells.num+1:end);
        
        CrNew_Na = GC.chemical.ion.Na(1:GC.Matrix.cells.num);
        CfNew_Na = GC.chemical.ion.Na(GC.Matrix.cells.num+1:end);       
        
        CrNew_Al = GC.chemical.ion.Al(1:GC.Matrix.cells.num);
        CfNew_Al = GC.chemical.ion.Al(GC.Matrix.cells.num+1:end);
        
        CrNew_Mg = GC.chemical.ion.Mg(1:GC.Matrix.cells.num);
        CfNew_Mg = GC.chemical.ion.Mg(GC.Matrix.cells.num+1:end);
        
        CrNew_H = GC.chemical.ion.H(1:GC.Matrix.cells.num);
        CfNew_H = GC.chemical.ion.H(GC.Matrix.cells.num+1:end);

        pHr = GC.chemical.pH(1:GC.Matrix.cells.num);
        pHf = GC.chemical.pH(GC.Matrix.cells.num+1:end);
        
        resp = Pr-pNew;
        resT = Tr-TNew;
        resNorm = max( norm(resp),norm(resT)); 
        
        nit = nit + 1;
        fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
    end

    [Gmfr,T] = getNewPropsP(Gmfr,GC,F,fracture); % couple C-H 
    [GT,TT] = getNewPropsT(GT,GC,F,fracture,cond_fluid,cond_solid); % couple C-T
    poroC = GC.rock.poro(:);
    permC = GC.rock.perm(:);
    quartzC = GC.chemical.quartz(:);
    ksparC = GC.chemical.kspar(:);
    albiteC = GC.chemical.albite(:); 
    phlogopiteC = GC.chemical.phlogopite(:);
    
    
    
    % time steps determine in terms of CFL number
    deltaL = physdim(1)/celldim(1)*0.6;
    dt =  deltaL/max(max(velocity.m),max(velocity.fr));
    dtc = deltaL/max(max(velocity.m),max(velocity.fr));
    iterc    = round(dtc/dt);
    
    [u_x,u_y,nnx,nny,node] = xfvm_2d(Gmfr, L, D, nx, ny, prNew, pfNew, TrNew, TfNew, Tri, Tfi, betaT, fl, F);
    
    if nit > maxits
        error('Newton solves did not converge')
    else % store solution
        sol(step+1)  = struct('time', convertTo(t, year), 'prNew', prNew,'pfNew', pfNew, 'pinj', pinj, 'Fluxinj', Fluxinj, 'Fluxout', Fluxout, 'TrNew',TrNew,'TfNew',TfNew,'Pwinj',Pwinj,'Pwout',Pwout,'Tout',Tout,'CrNew_SiO2',CrNew_SiO2,'CfNew_SiO2',CfNew_SiO2,'CrNew_K',CrNew_K,'CfNew_K',CfNew_K,'CrNew_Na',CrNew_Na,'CfNew_Na',CfNew_Na,'CrNew_Al',CrNew_Al,'CfNew_Al',CfNew_Al,'CrNew_Mg',CrNew_Mg,'CfNew_Mg',CfNew_Mg,'CrNew_H',CrNew_H,'CfNew_H',CfNew_H,'pHr',pHr,'pHf',pHf,'quartzC',quartzC,'ksparC',ksparC,'albiteC',albiteC,'phlogopiteC',phlogopiteC,'poroC',poroC,'permC',permC,'u_x',u_x,'u_y',u_y);
        waitbar(t/totTime,hwb);
    end
    
    substates.pressure = [stateP.prNew;stateP.pfNew;];
    substates.temperature = [stateT.TrNew;stateT.TfNew;];
    substates.SiO2  = [CrNew_SiO2;CfNew_SiO2];
    substates.K = [CrNew_K;CfNew_K];
    substates.Na = [CrNew_Na;CfNew_Na];
    substates.Al = [CrNew_Al;CfNew_Al];
    substates.Mg  = [CrNew_Mg;CfNew_Mg];
    substates.H = [CrNew_H;CfNew_H];
    substates.pH = [pHr;pHf];
    substates.quartzC = quartzC;
    substates.ksparC = ksparC;
    substates.albiteC = albiteC;
    substates.phlogopiteC = phlogopiteC;
    states(step) = {substates};

    if mod(step,10) == 0
        fileName = [baseFileName,'.mat'];
        save(fileName);
    end
    
toc
end
close(hwb)

%% ## License
% Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.
% 
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% 
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see <http://www.gnu.org/licenses/>.

% This project builds upon the original work of TU Delft and SINTEF ICT by adding new features and modifications. 
% All new code is released under the same GNU General Public License, version 3, or any later version.
% 
% This software is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% **Note:** 
% Certain functions or sections of the code may be licensed under different terms. 
% Please refer to the specific license details provided in the comments or documentation 
% of the corresponding functions or files.
% 
% Other licenses: please see details in their corresponding function descriptions.


