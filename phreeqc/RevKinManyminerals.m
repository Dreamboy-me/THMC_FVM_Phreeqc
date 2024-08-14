function [qrea_Si, qrea_K, qrea_Na, qrea_Al, qrea_Mg, qrea_H, pH, Mquartz, Mkspar, Malbite, Mphlogopite, dphi, dphiquartz, dphikspar, dphialbite, dphiphlogopite] = RevKinManyminerals(G,  prNew,  pfNew, TrNew, TfNew, rhof, dtc)

% Parameters
nc = G.cells.num;
Temp = [TrNew;TfNew];
press = [prNew./1.0e6;pfNew./1e6]; % Mpa
rhoq = G.chemical.rho.rock; % density，g/m^3
if isfield(G.rock,'poro')
    pv = poreVolume(G,G.rock); % m^3
else
    pv = G.cells.volumes;
end
M_SiO2 = G.chemical.Molemass.quartz;  % g/mol
M_Kspar = G.chemical.Molemass.kspar;
M_Albite = G.chemical.Molemass.albite;
M_Phlogopite = G.chemical.Molemass.phlogopite;

Siconc0 = G.chemical.ion.SiO2(:)./1000; % mol/m^3-->mol/Lw
Kconc0 = G.chemical.ion.K(:)./1000;
Naconc0 = G.chemical.ion.Na(:)./1000;
Alconc0 = G.chemical.ion.Al(:)./1000;
Mgconc0 = G.chemical.ion.Mg(:)./1000;
Hconc0 =  G.chemical.ion.H(:)./1000; % mol/m^3-->mol/Lw
pHconc0 = -log10(G.chemical.ion.H(:)./1e3); % mol/L  

quartzmoles  = G.chemical.moles.quartz(:); % mol/L
quartzmoles0 = G.chemical.moles.quartz(:);
ksparmoles  = G.chemical.moles.kspar(:);
ksparmoles0 = G.chemical.moles.kspar(:);
albitemoles  = G.chemical.moles.albite(:);
albitemoles0 = G.chemical.moles.albite(:);
phlogopitemoles  = G.chemical.moles.phlogopite(:);
phlogopitemoles0 = G.chemical.moles.phlogopite(:);

Siconc = zeros(nc,1);
Kconc = zeros(nc,1);
Naconc = zeros(nc,1);
Alconc = zeros(nc,1);
Mgconc = zeros(nc,1);
Hconc = zeros(nc,1);

dc_quartz = zeros(nc,1);
dc_kspar = zeros(nc,1);
dc_albite = zeros(nc,1);
dc_phlogopite = zeros(nc,1);

parfor i = 1:nc
    % Define input string as cell array of strings and combine them to a multiline string using the sprintf command.
    temperature = Temp(i);   % ºC
    pressure = press(i).*10; % atm
    density = rhof(i)./1000; % kg/m^3w --> kg/Lw
    Si = Siconc0(i); % mol/Lw
    K = Kconc0(i);
    Na = Naconc0(i);
    Al = Alconc0(i);
    Mg = Mgconc0(i);
    pH0 = pHconc0(i);
    
    reaction_pressure = press(i).*10; % atm
    reaction_temperature = Temp(i);
    
    quartzm = quartzmoles(i); % mol/Lw
    quartzm0 = quartzmoles0(i);
    ksparm = ksparmoles(i);
    ksparm0 = ksparmoles0(i);
    albitem = albitemoles(i);
    albitem0 = albitemoles0(i);
    phlogopitem = phlogopitemoles(i);
    phlogopitem0 = phlogopitemoles0(i);
    
    OUTphreeqSTRING = phreeqcSolvermanyminerals(temperature, pressure, density, Si, K, Na, Al, Mg, pH0, reaction_pressure, reaction_temperature, quartzm, quartzm0, ksparm, ksparm0, albitem, albitem0, phlogopitem, phlogopitem0, dtc);
    
    % Sum conc. of H+ in solution
    out_PHREEQC_H  = cell2mat(OUTphreeqSTRING(end,14));  
    Hconc(i) = out_PHREEQC_H;  % mol/kgw
    % Sum conc. of SiO2 in solution after all reaction time
    out_PHREEQC_Si = cell2mat(OUTphreeqSTRING(end,9));   
    Siconc(i) = out_PHREEQC_Si;  % mol/kgw
    out_PHREEQC_K = cell2mat(OUTphreeqSTRING(end,10));
    Kconc(i) = out_PHREEQC_K;
    out_PHREEQC_Na = cell2mat(OUTphreeqSTRING(end,11));
    Naconc(i) = out_PHREEQC_Na;
    out_PHREEQC_Al = cell2mat(OUTphreeqSTRING(end,12));
    Alconc(i) = out_PHREEQC_Al;
    out_PHREEQC_Mg = cell2mat(OUTphreeqSTRING(end,13));
    Mgconc(i) = out_PHREEQC_Mg;
    % pH in solution
    out_PHREEQC_pH = cell2mat(OUTphreeqSTRING(end,7));
    pH(i) = out_PHREEQC_pH;
    
    % Delta moles of quartz during the reaction process
    out_PHREEQC_quartz = cell2mat(OUTphreeqSTRING(2:end,20));  
    dc_quartz(i) = sum(out_PHREEQC_quartz).*(-1); % mol/Lw 
    out_PHREEQC_kspar = cell2mat(OUTphreeqSTRING(2:end,22));
    dc_kspar(i) = sum(out_PHREEQC_kspar).*(-1);
    out_PHREEQC_albite = cell2mat(OUTphreeqSTRING(2:end,24));
    dc_albite(i) = sum(out_PHREEQC_albite).*(-1);
    out_PHREEQC_phlogopite = cell2mat(OUTphreeqSTRING(2:end,26));
    dc_phlogopite(i) = sum(out_PHREEQC_phlogopite).*(-1);
    
end
 
qrea_Si = (Siconc - Siconc0 *1000./rhof).*rhof./dtc .* pv;  % mol/s 
qrea_K  = (Kconc - Kconc0 *1000./rhof).*rhof./dtc .* pv;
qrea_Na = (Naconc - Naconc0 *1000./rhof).*rhof./dtc .* pv;
qrea_Al = (Alconc - Alconc0 *1000./rhof).*rhof./dtc .* pv;
qrea_Mg = (Mgconc - Mgconc0 *1000./rhof).*rhof./dtc .* pv;
qrea_H  = (Hconc - Hconc0 *1000./rhof).*rhof./dtc .* pv;   % mol/s

% Mass of minerals in the rock undergoing dissolution/precipitation
dmquartz = dc_quartz.* pv *1.0e3 * M_SiO2; % g 
Mquartz = quartzmoles(:) - dc_quartz(:); % mol/L
dmkspar = dc_kspar.* pv *1.0e3 * M_Kspar;
Mkspar = ksparmoles(:) - dc_kspar(:);
dmalbite = dc_albite.* pv *1.0e3 * M_Albite;
Malbite = albitemoles(:) - dc_albite(:);
dmphlogopite = dc_phlogopite.* pv *1.0e3 * M_Phlogopite;
Mphlogopite = phlogopitemoles(:) - dc_phlogopite(:);

% poro change and src calculation 
dm = dmquartz + dmkspar + dmalbite + dmphlogopite;
dv = dm./rhoq;  % m^3  
dphi = dv./G.cells.volumes; % Porosity variation

dvquartz = dmquartz./rhoq;
dphiquartz = dvquartz./G.cells.volumes;
dvkspar = dmkspar./rhoq;
dphikspar = dvkspar./G.cells.volumes;
dvalbite = dmalbite./rhoq;
dphialbite = dvalbite./G.cells.volumes;
dvphlogopite = dmphlogopite./rhoq;
dphiphlogopite = dvphlogopite./G.cells.volumes;


end
