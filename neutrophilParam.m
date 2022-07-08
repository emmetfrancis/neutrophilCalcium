function p = neutrophilParam(p,condition)

% if ~isfield(p,'Ac')
%     p.t = (-1000:1000)';
%     t = p.t;
%     Amin = 0;
%     Amax = 240;
%     t0 = 20;
%     p.Ac = Amin + (Amax-Amin) ./ (1 + exp(-p.t/t0));
%     p.dAcdt = (Amax-Amin)*exp(-t/t0)./(t0*(1 + exp(-t/t0)).^2);
%     p.rhoIgG = 1000;
%     p.startTime = -100;
% end

F = .096480; % Faraday's constant with relevant units, Coulombs per umol
N_A = 6.022e17; % Avogadros number in relevant units, molecules per umol
R = 8.314e-3; % gas constant in mJ/(K*umol)
T = 293; % temperature in Kelvin
ERVolFraction = 0.05; %1
cytosolVolFraction = 0.63;
rhoFcR = 1700;
D = 0.01; % receptor diffusion coefficient in um^2 per second, order of mag
Ka_RL = 30; % 2D receptor-ligand binding affinity
ITAMFraction = 0.06;
Dcell = 8.75;
microSA = 2.1;
vol_tot = (4/3)*pi*(Dcell/2)^3;
vol_i = cytosolVolFraction*vol_tot*1e-15; %exclude organelles
vol_ER = ERVolFraction*vol_tot*1e-15; %approx, can adjust
vol_iOther = 1e-12;
nu_Src = 0.5;%0.1;
K_Src = 10;
k_dePhos = 0.1;
alpha_PLC = 2.5e-4;%1e-4;%5e5*2.8e-5/(N_A*vol_i); %s^-1, Lemon et al 2003
K_PLC = 0.4; %uM, Lemon et al 2003
r_r = .015; %s^-1, Lemon et al 2003
k_degIP3 = 1.25; %s^-1, Fink et al 1999
PIP2_tot = 2e7;%5e7; % #molecules, Lemon et al 2003
P_IP3 = 2880; %uM/s, Fink et al 1999, check these units, key to osc
K_iIP3R = 0.1; %uM, Fink et al 1999
K_aIP3R = 0.17; %uM, Fink et al 1999
A_h = 1.4; % 1/(uM s), Bennett et al 2005, sets timescale of channel deactivation
K_dh = 0.12; % deactivation of channels due to calcium

% nu_leakER = 0.1; % could be ~0.3 to 1.5 s^-1, refer to Camello et al 2002
nu_SERCA = 51.8242;%10e-12 / (2*F*vol_iOther); %uM/s, converted from picoAmperes, Yang et al 2003
K_SERCA = 0.5;%1; %uM, Miller and Carsten 1997, key to deterimining if osc occur
h0 = K_dh/(.08 + K_dh);
IP30 = 0;
c_i0 = .1;
% nu_leakER = nu_SERCA*c_i0^2/((K_SERCA^2 + c_i0^2)*(200-c_i0)) - ...
%     P_IP3*IP30.^3 .* c_i0.^3 .* h0.^3 ./ ((IP30 + K_iIP3R).^3 .* (c_i0 + K_aIP3R).^3.*200)
nu_leakER = 0.0099712;
g_SOCE = (1e-12)/200; %in units of Amperes per mV
c_ref = 70;%0.1; % uM, from measurements of soce in immune cells
j_0PMCA = 0.8387385;%0.025*5e-12 / (2*F*vol_iOther); %uM/s, converted from picoAmperes, estimated in Kapela et al
K_mPMCA = 0.17; %uM, O'Donnell and Owen 1994

% fast buffering parameters
B_itot = 300;%760; 
K_iBuffer = 0.5;
B_iBAPTA = .1e3;
K_iBAPTA = 0.7;
B_ERtot = 60e3; % CSQN from Kapela et al
K_ERBuffer = 500;
% Voltage parameters
V_PMN = -60; % membrane voltage in mV of human neutrophil

% other parameters...
g_leakPM = 0;
dtSignal = 10;

Amax = p.Amax;
tShift = 100;
t0 = 20;

pVec = zeros(35,1);
pVec(1) = Dcell;
pVec(2) = microSA;
pVec(3) = cytosolVolFraction;
pVec(4) = ERVolFraction;
pVec(5) = p.rhoIgG;
pVec(6) = rhoFcR;
pVec(7) = ITAMFraction;
pVec(8) = D;
pVec(9) = dtSignal;
pVec(10) = Ka_RL;
pVec(11) = nu_Src;
pVec(12) = K_Src;
pVec(13) = k_dePhos;
pVec(14) = alpha_PLC;
pVec(15) = K_PLC;
pVec(16) = r_r;
pVec(17) = k_degIP3;
pVec(18) = PIP2_tot;
pVec(19) = P_IP3;
pVec(20) = K_iIP3R;
pVec(21) = K_aIP3R;
pVec(22) = K_dh;
pVec(23) = A_h;
pVec(24) = nu_SERCA;
pVec(25) = K_SERCA;
pVec(26) = nu_leakER;
pVec(27) = g_SOCE;
pVec(28) = c_ref;
pVec(29) = j_0PMCA;
pVec(30) = K_mPMCA;
pVec(31) = B_itot;
pVec(32) = K_iBuffer;
pVec(33) = B_ERtot;
pVec(34) = K_ERBuffer;
pVec(35) = V_PMN; 
pVec(36) = Amax;
pVec(37) = tShift;
pVec(38) = t0;


p.pVec = pVec;

p.IP30 = IP30;

% flags
switch condition
    case 'Standard'
        p.c_EC = 1260; % 1.26 mM extracellular calcium
        p.flag_pITAM = 1;
        p.flag_IP3Prod = 1;
        p.flag_IP3R = 1;
        p.flag_otherR = 0;
        p.flag_SERCA = 1;
        p.flag_SOCE = 1;
        p.flag_PMCA = 1;
        p.flag_BAPTA = 0;
    case 'Ca Free'
        p.c_EC = 0;%1e-6; 
        p.flag_pITAM = 1;
        p.flag_IP3Prod = 1;
        p.flag_IP3R = 1;
        p.flag_otherR = 0;
        p.flag_SERCA = 1;
        p.flag_SOCE = 1;
        p.flag_PMCA = 1;
        p.flag_BAPTA = 0;
    case 'Thaps'
        p.c_EC = 0;%1e-6; 
        p.flag_pITAM = 1;
        p.flag_IP3Prod = 1;
        p.flag_IP3R = 1;
        p.flag_otherR = 0;
        p.flag_SERCA = 0.1;
        p.flag_SOCE = 1;
        p.flag_PMCA = 1;
        p.flag_BAPTA = 0;
    case 'BAPTA'
        p.c_EC = 1260; % 1.26 mM extracellular calcium
        p.flag_pITAM = 1;
        p.flag_IP3Prod = 1;
        p.flag_IP3R = 1;
        p.flag_otherR = 0;
        p.flag_SERCA = 1;
        p.flag_SOCE = 1;
        p.flag_PMCA = 1;
        p.flag_BAPTA = 1;
    case 'PLCInhib'
        p.c_EC = 1260; % 1.26 mM extracellular calcium
        p.flag_pITAM = 1;
        p.flag_IP3Prod = 1;
        p.flag_IP3R = 1;
        p.flag_otherR = 0;
        p.flag_SERCA = 1;
        p.flag_SOCE = 1;
        p.flag_PMCA = 1;
        p.flag_BAPTA = 0;
        p.pVec(10) = p.pVec(10) * 0.01;
end