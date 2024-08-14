function G = incompTPFACWell(state_C, G, F, fracture, fluid, velocity, prNew,  pfNew, TrNew, TfNew, rhof, mu, Diff, dtc, Cinj, bcC, W) 


[qrea_Si, qrea_K, qrea_Na, qrea_Al, qrea_Mg, qrea_H, pH, Mquartz, Mkspar, Malbite, Mphlogopite, dphi, dphiquartz, dphikspar, dphialbite, dphiphlogopite] = RevKinManyminerals(G,  prNew,  pfNew, TrNew, TfNew, rhof, dtc); 

% Mineral content changes
G.chemical.quartz = G.chemical.quartz(:) - dphiquartz(:)./(1-G.rock.poro(:)); % Variation in the volume fraction
G.chemical.kspar = G.chemical.kspar(:)   - dphikspar(:)./(1-G.rock.poro(:));
G.chemical.albite = G.chemical.albite(:) - dphialbite(:)./(1-G.rock.poro(:));
G.chemical.phlogopite = G.chemical.phlogopite(:) - dphiphlogopite(:)./(1-G.rock.poro(:));
G.rock.poro = G.rock.poro(:) + dphi(:);  

G.chemical.pH  = pH(:);
G.chemical.moles.quartz  = Mquartz; % mol/L  
G.chemical.moles.kspar   = Mkspar;
G.chemical.moles.albite  = Malbite;
G.chemical.moles.phlogopite  = Mphlogopite;

% for effective diffusion & TI update 
[G,T] = getNewPropsC(G, F, fracture, Diff); 
 
src_Si = addSource([], G.cells.indexMap, qrea_Si); %  mol/s
src_K  = addSource([], G.cells.indexMap, qrea_K);
src_Na = addSource([], G.cells.indexMap, qrea_Na);
src_Al = addSource([], G.cells.indexMap, qrea_Al);
src_Mg = addSource([], G.cells.indexMap, qrea_Mg);
src_H  = addSource([], G.cells.indexMap, qrea_H); %  mol/s 

CrOld_SiO2   =  G.chemical.ion.SiO2(1:G.Matrix.cells.num); % mol/m^3
CfOld_SiO2   =  G.chemical.ion.SiO2(G.Matrix.cells.num+1:end);
stateC = incompTPFA_C(state_C, G, T, fluid, velocity, prNew, CrOld_SiO2, CfOld_SiO2, rhof, mu, dtc, Cinj, 'bc', bcC, 'src', src_Si, 'Wells', W, 'use_trans', true);
CrNew_SiO2 = stateC.CrNew;
CfNew_SiO2 = stateC.CfNew;
G.chemical.ion.SiO2 = [CrNew_SiO2;CfNew_SiO2]; % mol/m^3

CrOld_K   =  G.chemical.ion.K(1:G.Matrix.cells.num); % mol/m^3
CfOld_K   =  G.chemical.ion.K(G.Matrix.cells.num+1:end);
stateC = incompTPFA_C(state_C, G, T, fluid, velocity, prNew, CrOld_K, CfOld_K, rhof, mu, dtc, Cinj, 'bc', bcC, 'src', src_K, 'Wells', W, 'use_trans', true);
CrNew_K = stateC.CrNew;
CfNew_K = stateC.CfNew;
G.chemical.ion.K = [CrNew_K;CfNew_K]; % mol/m^3

CrOld_Na   =  G.chemical.ion.Na(1:G.Matrix.cells.num); % mol/m^3
CfOld_Na   =  G.chemical.ion.Na(G.Matrix.cells.num+1:end);
stateC = incompTPFA_C(state_C, G, T, fluid, velocity, prNew, CrOld_Na, CfOld_Na, rhof, mu, dtc, Cinj, 'bc', bcC, 'src', src_Na, 'Wells', W, 'use_trans', true);
CrNew_Na = stateC.CrNew;
CfNew_Na = stateC.CfNew;
G.chemical.ion.Na = [CrNew_Na;CfNew_Na]; % mol/m^3

CrOld_Al   =  G.chemical.ion.Al(1:G.Matrix.cells.num); % mol/m^3
CfOld_Al   =  G.chemical.ion.Al(G.Matrix.cells.num+1:end);
stateC = incompTPFA_C(state_C, G, T, fluid, velocity, prNew, CrOld_Al, CfOld_Al, rhof, mu, dtc, Cinj, 'bc', bcC, 'src', src_Al, 'Wells', W, 'use_trans', true);
CrNew_Al = stateC.CrNew;
CfNew_Al = stateC.CfNew;
G.chemical.ion.Al = [CrNew_Al;CfNew_Al]; % mol/m^3

CrOld_Mg   =  G.chemical.ion.Mg(1:G.Matrix.cells.num); % mol/m^3
CfOld_Mg   =  G.chemical.ion.Mg(G.Matrix.cells.num+1:end);
stateC = incompTPFA_C(state_C, G, T, fluid, velocity, prNew, CrOld_Mg, CfOld_Mg, rhof, mu, dtc, Cinj, 'bc', bcC, 'src', src_Mg, 'Wells', W, 'use_trans', true);
CrNew_Mg = stateC.CrNew;
CfNew_Mg = stateC.CfNew;
G.chemical.ion.Mg = [CrNew_Mg;CfNew_Mg]; % mol/m^3

CrOld_H   =  G.chemical.ion.H(1:G.Matrix.cells.num); % mol/m^3
CfOld_H   =  G.chemical.ion.H(G.Matrix.cells.num+1:end);
stateC = incompTPFA_C(state_C, G, T, fluid, velocity, prNew, CrOld_H, CfOld_H, rhof, mu, dtc, Cinj, 'bc', bcC, 'src', src_H, 'Wells', W, 'use_trans', true);
CrNew_H = stateC.CrNew; 
CfNew_H = stateC.CfNew;
G.chemical.ion.H = [CrNew_H;CfNew_H]; % mol/m^3

end