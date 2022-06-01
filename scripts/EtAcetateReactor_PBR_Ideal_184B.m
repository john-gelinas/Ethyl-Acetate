function reaction = EtAcetateReactor_PBR_Ideal_184B(V,N,T,P)

%T is in Kelvin
%P is in atm
R = 1.987;%[cal/mol*K]
packingfraction = 0.5*6.31; %g/cm^3 (packing fraction*CuO density)
%set values for molar flowarates at current step
na = N(1);
nr = N(2);
ns = N(3);
nea = N(4);
ntot = sum(N);
%calc partial pressures at given step
pa = P*na/ntot;%[atm]
pea = P*nea/ntot;%[atm]
ps = P*ns/ntot;%[atm]
%calc equilibrium values
K = exp((-3127/(T))+5.32);%[1/atm]
k = exp((-16310/(R*T))+10.95);%[mol/hr*g-cat]
KA = exp((5890/(R*T))-6.40);%[1/atm]
Kea = exp((11070/(R*T))-9.40);%[1/atm]
KS = exp((6850/(R*T))-7.18);%[1/atm]
%calc reaction rate for ethanol on catalyst
ra = k*KA*((pa^2) - (pea*(ps^2)/K))/((1 + KA*pa + Kea*pea + KS*ps)^2);%[mol/hr*g-cat]
ra = ra*packingfraction; %[mol/hr*cm^3]
ra = ra*1000; %[mol/hr*L]
dH2 = 2*ra;
% dAcetaldehyde = ra;
dEtAcetate = ra;
dEther = (1/10)*ra;
dH2O = dEther;
dAcetaldehyde = 0;
dEthanol = -2*dEtAcetate-2*dEther;

reaction = [dEthanol;dAcetaldehyde;dH2;dEtAcetate;dEther;dH2O];
end