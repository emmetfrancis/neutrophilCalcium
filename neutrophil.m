function [dy,fluxes] = neutrophil(t,y,p)
% parameter vector p - minimally contains contact area curve info and IgG
% density
% p{1} = tAc (time vector for Ac curve)
% p{2} = Ac (contact area values assoc with time vector tAc)
% p{3} = dAcdt (derivative of contact area over time)
% p{4} = rhoIgG (IgG density on target)
% p{5} = rhoFcR (receptor density on cell)
% opt. p{6} = sigmaParam (for population tests)

% state variables (concentrations in micromolar)
% y(1): c_i (free cytosolic calcium concentration)
% y(2): c_ER (free endoplasmic reticulum calcium concentration)
% y(3): R_bound (number of bound receptors in contact region)
% y(4): pITAM (number of phosphorylated ITAM domains in the
% contact region
% y(5): PIP2 (number of PIP2 in contact region)
% y(6): IP3 (concentration of IP3 in cytosol)
% y(7): S (concentration of additional second messenger in cytosol)
% y(8): h (fraction of IP3R channels not deactivated by calcium)

% define variables
c_i = y(1);
c_ER = y(2);
R_bound = y(3);
pITAM = y(4);
PIP2 = y(5);
IP3 = y(6);
S = y(7);
h = y(8);

dy = zeros(size(y));
    
% define parameters
F = .096480; % Faraday's constant with relevant units, Coulombs per umol
N_A = 6.022e17; % Avogadros number in relevant units, molecules per umol
R = 8.314e-3; % gas constant in mJ/(K*umol)
T = 293; % temperature in Kelvin
tAc = p.t;
Ac = p.Ac;
dAcdt = p.dAcdt;
% rhoIgG = p.rhoIgG;
c_EC = p.c_EC;

% pVec contains all the independent parameters for population tests
pVec = p.pVec;
Dcell = pVec(1);
microSA = pVec(2);
cytosolVolFraction = pVec(3);
ERVolFraction = pVec(4);
rhoIgG = pVec(5);
rhoFcR = pVec(6);
ITAMFraction = pVec(7);
D = pVec(8);
dtSignal = pVec(9);
Ka_RL = pVec(10);
nu_Src = pVec(11);
K_Src = pVec(12);
k_dePhos = pVec(13);
alpha_PLC = pVec(14);
K_PLC = pVec(15);
r_r = pVec(16);
k_degIP3 = pVec(17);
PIP2_tot = pVec(18);
P_IP3 = pVec(19);
K_iIP3R = pVec(20);
K_aIP3R = pVec(21);
K_dh = pVec(22);
A_h = pVec(23);
nu_SERCA = pVec(24);
K_SERCA = pVec(25);
nu_leakER = pVec(26);
g_SOCE = pVec(27);
c_ref = pVec(28);
j_0PMCA = pVec(29);
K_mPMCA = pVec(30);
B_itot = pVec(31);
K_iBuffer = pVec(32);
B_ERtot = pVec(33);
K_ERBuffer = pVec(34);
V_PMN = pVec(35); 

% calculate other values assoc with parameters
Atot = 4*pi*(Dcell/2)^2 * microSA; %microscopic membrane area including folds
Rtot = Atot*rhoFcR;
dL = 1; % characteristic length for receptor diffusion
Dbar = sqrt(4*pi)*D/dL;
vol_tot = (4/3)*pi*(Dcell/2)^3;
vol_i = cytosolVolFraction*vol_tot*1e-15; %exclude organelles
vol_ER = ERVolFraction*vol_tot*1e-15; %approx, can adjust
IP30 = 0;

B_iBAPTA = 1e3;
K_iBAPTA = 0.7;
E_Ca = (R*T/(2*F)) * log(c_EC/c_i);
if E_Ca < V_PMN
    E_Ca = V_PMN;
end
%
% E_Ca0 = (R*T/(2*F)) * log(1260/c_i);
% c_EC0 = 1260;
% c0Adj = .08;% - .05;
% g_SOCE = j_0PMCA*(c0Adj/(c0Adj+K_mPMCA))*(2*F*vol_i*(1+200/c_ref)/(E_Ca0-V_PMN));

% % NCX parameters
% IbarNCX = 0.1*(4.5+1); % [A/F] (Original 4.5 Grandi-Bers) 0.5
% g_NCX = (IbarNCX*10e-12) / (2*F*vol_i); % converted to [uM/s]
% Qcorr = (T-310)/10;
% Nai = 10; %[mM], need to determine amt for sure in pmn
% Nao = 142.4; %[mM] in HBSS+
% V = V_PMN;
% RT_F = R*T/F;
% Kmc_i = 6.59e-3; % [mM] 3.59e-3; original
% Kmc_EC = 1.3; % [mM]
% KmNai = 12.29; % [mM]
% KmNao = 87.5; % [mM]
% ksat = 0.32; % [none]
% nu = 0.27; % [none]
% Kdact = 0.050e-3; % [mM]0.150e-3 original
% Q10NCX = 1.57; % [none]
% Ka = 1/(1+(Kdact/c_i0)^2);
% s1 = exp(nu*V/RT_F)*Nai^3*c_EC;
% s2 = exp((nu-1)*V/RT_F)*Nao^3*c_i0;
% s3 = Kmc_i*Nao^3*(1+(Nai/KmNai)^3) + KmNao^3*c_i0*(1+c_i0/Kmc_i)+Kmc_EC*Nai^3+Nai^3*c_EC+Nao^3*c_i0;
% j0_NCX = g_NCX*Q10NCX^Qcorr*Ka*(s1-s2)/s3/(1+ksat*exp((nu-1)*V/RT_F));

% g_leakPM = (1/(c_EC0-.08)).*(j_0PMCA*(c0Adj/(c0Adj+K_mPMCA)) - j0_NCX - g_SOCE*(E_Ca0 - V_PMN)./(2*F*vol_i*(1+200/c_ref)));
% g_leakPM = (1/(c_EC0-.08)).*(j_0PMCA*(c0Adj/(c0Adj+K_mPMCA)) - g_SOCE*(E_Ca0 - V_PMN)./(2*F*vol_i*(1+200/c_ref)));

% load flags
flag_pITAM = p.flag_pITAM;
flag_IP3Prod = p.flag_IP3Prod;
flag_IP3R = p.flag_IP3R;
flag_SERCA = p.flag_SERCA;
flag_SOCE = p.flag_SOCE;
flag_PMCA = p.flag_PMCA;
flag_BAPTA = p.flag_BAPTA;
% flag_NCX = 0;%1;%.01;

%% start calculations
%patch spread over past dtSignal s
A_front = interp1(tAc,Ac,t) - interp1(tAc,Ac,t-dtSignal);

% receptor system
Ac = interp1(tAc,Ac,t);
dAcdt = interp1(tAc,dAcdt,t);

R_bound = R_bound./Ac;
R_cup = Ka_RL*R_bound./(rhoIgG-R_bound) + R_bound;
R_free = R_cup - R_bound;
F_R = 1./(1 + rhoIgG*Ka_RL./(Ka_RL + R_free).^2);
firstTerm = -R_cup.*((Ac./(Atot-Ac)).*(Dbar*sqrt(Ac) + dAcdt) + dAcdt);
secondTerm = -Dbar*sqrt(Ac).*(R_free - Rtot./(Atot-Ac)) + Rtot.*dAcdt./(Atot-Ac);
dR = (F_R./Ac) .* (firstTerm + secondTerm);
% convert to change in bound number over time
dy3Orig = Ac.*dR.*(1./F_R - 1) + R_bound.*dAcdt;
% alt definition
N_bound = R_bound*Ac;
N_tot = Rtot;
N_cup = N_bound*(1 + Ka_RL*Ac./(rhoIgG*Ac - N_bound));
F_RAlt = 1./(1 + Ka_RL*rhoIgG*Ac.^2./((rhoIgG*Ac-N_bound).^2));
rho_RBody = (N_tot-N_cup)/(Atot-Ac);
diffusionTerm = D*2*pi*sqrt(Ac/pi)*(rho_RBody*Ac - (N_cup-N_bound))./(dL*Ac); 
otherTerms = dAcdt*rho_RBody - dAcdt*Ka_RL./((rhoIgG*Ac-N_bound).^2);
dy(3) = F_RAlt*(diffusionTerm + otherTerms);
if isnan(dy(3))
    dy(3) = 0;
end

%change in pITAM over time
totITAM = ITAMFraction*R_bound.*Ac;
ITAM = totITAM - pITAM;
nu_SrcInner = 0;
dy(4) = flag_pITAM*((A_front*nu_Src*ITAM./Ac +...
    (Ac-A_front).*nu_SrcInner.*ITAM./Ac).*ITAM./(K_Src*Ac + ITAM) - k_dePhos*pITAM);

%PIP2 <-> IP3
% IP30 = 3;
dy(5) = -(alpha_PLC*pITAM.* (c_i/(K_PLC+c_i)) + r_r) .* PIP2 -...
    r_r*vol_i*N_A*IP3 + r_r*PIP2_tot; %PIP2 eq
dy(6) = flag_IP3Prod*(alpha_PLC*(pITAM).* (c_i/(K_PLC+c_i)).* (PIP2/(vol_i*N_A)) - k_degIP3*(IP3-IP30)); %IP3 eq

%ER calcium fluxes
j_IP3R = flag_IP3R*(P_IP3*IP3.^3 .* c_i.^3 .* h.^3 ./ ((IP3 + K_iIP3R).^3 .* (c_i + K_aIP3R).^3)) .*...
    (c_ER - c_i)./200;%c_ER; %IP3R flux in uM/s
dy(8) = A_h*(K_dh - (c_i+K_dh).*h); % channel deactivation
j_leakER = nu_leakER * (c_ER - c_i);
if flag_SERCA == 0
    j_SERCA = (.1+.9*exp(-t/100))*nu_SERCA*c_i.^2./(K_SERCA.^2 + c_i.^2);
else
    j_SERCA = flag_SERCA*nu_SERCA*c_i.^2./(K_SERCA.^2 + c_i.^2);
end

%PM calcium fluxes
% j_SOCE = flag_SOCE*g_SOCE*1e-12*(E_Ca - V_PMN)./(2*F*vol_i*(1+c_ER/c_ref));
j_SOCE = flag_SOCE*g_SOCE*(E_Ca - V_PMN)./(2*F*vol_i*(1+(c_ER/c_ref)^4));
if j_SOCE < 0
    j_SOCE = 0;
end
cAdj = c_i;% - .05;
j_PMCA = flag_PMCA*j_0PMCA*cAdj./(cAdj + K_mPMCA);
if j_PMCA < 0
    j_PMCA = 0;
end
% j_leakPM = flag_PMCA*g_leakPM.*(c_EC - c_i);
j_leakPM = 0;

% % NCX
% Ka = 1/(1+(Kdact/c_i)^2);
% s1 = exp(nu*V/RT_F)*Nai^3*c_EC;
% s2 = exp((nu-1)*V/RT_F)*Nao^3*c_i;
% s3 = Kmc_i*Nao^3*(1+(Nai/KmNai)^3) + KmNao^3*c_i*(1+c_i/Kmc_i)+Kmc_EC*Nai^3+Nai^3*c_EC+Nao^3*c_i;
% j_NCX = flag_NCX*g_NCX*Q10NCX^Qcorr*Ka*(s1-s2)/s3/(1+ksat*exp((nu-1)*V/RT_F));

%calcium balances
f_i = 1/(1 + B_itot*K_iBuffer./((K_iBuffer+c_i).^2) + ...
    flag_BAPTA*B_iBAPTA*K_iBAPTA./((K_iBAPTA+c_i).^2));
f_ER = 1/(1 + B_ERtot*K_ERBuffer./((K_ERBuffer+c_ER).^2));
% dy(1) = f_i*(j_IP3R + j_other + j_leakER - j_SERCA + j_SOCE + j_leakPM - j_PMCA + j_NCX);
j_other = 0;
dy(1) = f_i*(j_IP3R + j_other + j_leakER - j_SERCA + j_SOCE + j_leakPM - j_PMCA);


dy(2) = f_ER*(vol_i/vol_ER)*(j_SERCA - j_IP3R - j_other - j_leakER);

% fluxes = [j_IP3R, j_other, j_leakER, j_SERCA, j_SOCE, j_leakPM, j_PMCA, j_NCX];
fluxes = [j_IP3R, j_other, j_leakER, j_SERCA, j_SOCE, j_leakPM, j_PMCA];

end