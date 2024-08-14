function G = chemical_initial(G, quartz, kspar, albite, phlogopite, clinochlore, varargin)
% ****** Volume fraction in the mineral phase ******
% multi components
nc = G.cells.num;
G.chemical.quartz  = quartz.*ones(nc,1); 
G.chemical.kspar   = kspar.*ones(nc,1);
G.chemical.albite  = albite.*ones(nc,1);
G.chemical.phlogopite  = phlogopite.*ones(nc,1);
G.chemical.clinochlore  = clinochlore.*ones(nc,1); 
G.chemical.unresponsive = (1 - quartz - kspar - albite - phlogopite).*ones(nc,1);

% ****** Density of rock and minerals ******
G.chemical.rho.rock    = 2.65e6; % g/m^3
G.chemical.rho.quartz  = 2.65e6; % g/m^3
G.chemical.rho.kspar   = 2.55e6;
G.chemical.rho.albite  = 2.62e6;
G.chemical.rho.phlogopite  = 2.78e6;
G.chemical.rho.clinochlore = 2.65e6;

% ****** Concentration under the equilibrium state of mineral ions ******
G.chemical.ion.SiO2  = 6.28e-3.*ones(nc,1)*1e3; % mol/kgw 
G.chemical.ion.K  = 2.89e-5.*ones(nc,1)*1e3;
G.chemical.ion.Na = 4.21e-4.*ones(nc,1)*1e3; 
G.chemical.ion.Al = 4.50e-4.*ones(nc,1)*1e3;  
G.chemical.ion.Mg = 4.86e-5.*ones(nc,1)*1e3;  
G.chemical.ion.H  = 1.0e-7.*ones(nc,1)*1e3;

% ################################ 
% Assumption: soil/rock is 15% Kspar, 32.1% Albite, 34.1% Quarz, 16.8% muscovite, and 2% Clinochlore in 1 mm spheres (radius 0.5 mm)
% Assumption: density of rock and (Kspar, Albite, Quarz, muscovite, and Clinochlore ) is 2650 kg/m^3 = 2.65 kg/L; 
% Assumption: pore space [-] =  XXX     # personal understanding--> which is the volume of solid skeleton in the rock of the geothermal
% GFW Kspar 0.278 kg/mol; GFW Albite 0.262 kg/mol; GFW Quartz 0.06 kg/mol; GFW muscovite 0.398 kg/mol; GFW Clinochlore 0.554 kg/mol;

% ****** The mass of rock and different types of minerals per liter of pore space ******
rockmass_per_litre = (1-G.rock.poro).* 2.65./G.rock.poro; % kg/L,  rockmass per litre pore space
quartzmass_per_litre = rockmass_per_litre .* G.chemical.quartz; % kg/L,  Quartzmass per litre  pore space
ksparmass_per_litre  = rockmass_per_litre .* G.chemical.kspar;
albitemass_per_litre = rockmass_per_litre .* G.chemical.albite;
phlogopitemass_per_litre = rockmass_per_litre .* G.chemical.phlogopite;
clinochloremass_per_litre = rockmass_per_litre .* G.chemical.clinochlore;

% ****** The molar amount of different types of minerals per liter of pore space, i.e., m/m0 ******
G.chemical.Molemass.quartz = 60.08;  
G.chemical.Molemass.kspar  = 278.33; 
G.chemical.Molemass.albite = 262.22; 
G.chemical.Molemass.phlogopite = 417.26; 
G.chemical.Molemass.clinochlore = 555.80; 

% Quartz: moles per liter pore space
quartzmoles_per_litre = quartzmass_per_litre ./(G.chemical.Molemass.quartz*1.0e-3); % mol Quartz/kgw 
G.chemical.moles.quartz = quartzmoles_per_litre(:); % mol/L

% Kspar: moles per liter pore space
ksparmoles_per_litre = ksparmass_per_litre ./(G.chemical.Molemass.kspar*1.0e-3);
G.chemical.moles.kspar = ksparmoles_per_litre(:); % mol/L

% Albite: moles per liter pore space
albitemoles_per_litre = albitemass_per_litre ./(G.chemical.Molemass.albite*1.0e-3);
G.chemical.moles.albite = albitemoles_per_litre(:); % mol/L

% Phlogopite: moles per liter pore space
phlogopitemoles_per_litre = phlogopitemass_per_litre ./(G.chemical.Molemass.phlogopite*1.0e-3);
G.chemical.moles.phlogopite = phlogopitemoles_per_litre(:); % mol/L

% Clinochlore: moles per liter pore space
clinochloremoles_per_litre = clinochloremass_per_litre ./(G.chemical.Molemass.clinochlore*1.0e-3);
G.chemical.moles.phlogopite = clinochloremoles_per_litre(:); % mol/L

end