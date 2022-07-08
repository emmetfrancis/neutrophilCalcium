function y0 = neutrophilInit(p)
    
    % load from struct
    Ac = p.Ac;
    t = p.t;
    startTime = p.startTime;
    pVec = p.pVec;
    rhoIgG = pVec(5);%p.rhoIgG;
    rhoFcR = pVec(6);
    Ka = pVec(10);
    PIP2_tot = pVec(18);

    % calculate y0
    [~,zeroIdx] = min(abs(t-startTime));
    bVal = Ka+rhoIgG-rhoFcR;
    R_free0 = -.5*(bVal)+.5*sqrt(bVal.^2+4*rhoFcR*Ka); %unbound receptor density in initial cup
    R_bound0 = rhoFcR - R_free0;
    fun = @(x) root3d(x,p);
    x = fsolve(fun,[.08,200,0.6],p);
    c_i0 = x(1);
    c_ER0 = x(2);
    h0 = x(3);
    
    y0 = [c_i0; c_ER0; R_bound0*Ac(zeroIdx); 0; PIP2_tot; 0; 0; h0];
    
end

function F = root3d(x,p)

pVec = p.pVec;
Dcell = pVec(1);
cytosolVolFraction = pVec(3);
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
V_PMN = pVec(35);
IP30 = p.IP30;
c_EC = p.c_EC;
if c_EC<.1
    c_EC = 1260;
end
N_A = 6.022e17; % Avogadros number in relevant units, molecules per umol
R = 8.314e-3; % gas constant in mJ/(K*umol)
T = 293; % temperature in Kelvin
F = .096480; % Faraday's constant with relevant units, Coulombs per umol
vol_tot = (4/3)*pi*(Dcell/2)^3;
vol_i = cytosolVolFraction*vol_tot*1e-15;
flag_IP3R = p.flag_IP3R;
flag_SERCA = 1;%p.flag_SERCA;
flag_SOCE = p.flag_SOCE;
flag_PMCA = p.flag_PMCA;

%ER calcium fluxes
j_IP3R = flag_IP3R*(P_IP3*IP30.^3 .* x(1).^3 .* x(3).^3 ./ ((IP30 + K_iIP3R).^3 .* (x(1) + K_aIP3R).^3)) .*...
    (x(2) - x(1))./200;%c_ER; %IP3R flux in uM/s
j_leakER = nu_leakER * (x(2) - x(1));
j_SERCA = flag_SERCA*nu_SERCA*x(1).^2./(K_SERCA.^2 + x(1).^2);

%PM calcium fluxes
E_Ca = (R*T/(2*F)) * log(c_EC/x(1));
if E_Ca < V_PMN
    E_Ca = V_PMN;
end
j_SOCE = flag_SOCE*g_SOCE*(E_Ca - V_PMN)./(2*F*vol_i*(1+(x(2)/c_ref)^4));
if j_SOCE < 0
    j_SOCE = 0;
end
cAdj = x(1);% - .07;
j_PMCA = flag_PMCA*j_0PMCA*cAdj./(cAdj + K_mPMCA);
% j_leakPM = flag_PMCA*g_leakPM.*(c_EC - x(1));

% find SS values
F(1) = j_SERCA - j_IP3R - j_leakER;
F(2) = j_SOCE - j_PMCA;
F(3) = K_dh - (x(1)+K_dh).*x(3);

end