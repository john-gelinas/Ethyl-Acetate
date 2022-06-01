function NPV = Sensitivity_Hysys_EthylAcetate(Price,EnterpriseRate,InterestRate,IncomeTaxRate,FCImultiplier,SUmultiplier,WCmultiplier)
% %Base Values%%
% Price.Eth = 500;%[$/MT]
% Price.DE = 800;%[$/MT]
% Price.EA = 1100;%[$/MT]
% Price.Fuel = 2.5;%[$/MM BTU]
% Price.CO2 = 40; %$/MT
% Price.Cat = 15000; %$/MT
% EnterpriseRate = 0.1;
% InterestRate = 0.05;
% IncomeTaxRate = 0.27;
% FCImultiplier = 1;
% SUmultiplier = 1;
% WCmultiplier = 1;
% Price.Electricity = 0.08; % $/kw hr

x = 0.45;
EAflow = 50000; %MT/yr
Ethflow = 59514; %MT/yr
DEflow = 3936.24; %MT/yr
Hflow = 4480.56; %MT/yr

MS = 1650;

%heat exchangers
HX1a = 1110; %ft^2
HX2a = 594; %ft^2
HX3a = 2023; %ft^2
HX4a = 1404; %ft^2
HX5a = 116; %ft^2
HX6a = 4; %ft^2
HX1Fc = (0.85+0)*3.75;
HX2Fc = (1.35+0)*3.75;
HX3Fc = (0.85+0)*3.75;
HX4Fc = (0.85+0)*3.75;
HX5Fc = (0.85+0)*3.75;
HX6Fc = (0.85+0)*3.75;

HX1cost = MS/280*101.3*HX1a^(0.65)*(2.29+HX1Fc); %$
HX2cost = MS/280*101.3*HX2a^(0.65)*(2.29+HX2Fc); %$
HX3cost = MS/280*101.3*HX3a^(0.65)*(2.29+HX3Fc); %$
HX4cost = MS/280*101.3*HX4a^(0.65)*(2.29+HX4Fc); %$
HX5cost = MS/280*101.3*HX5a^(0.65)*(2.29+HX5Fc); %$
HX6cost = MS/280*101.3*HX6a^(0.65)*(2.29+HX6Fc); %$

HX2duty = 3.16; %MW
HX3duty = 1.2; %MW
HX4duty = 4.07*3600*8400/1055; %MM BTU/year
HX5duty = 0.76*3600*8400/1055; %MM BTU/year
HX6duty = 0.01*3600*8400/1055; %MM BTU/year

HX2fuel = HX2duty*1.13; %MW
HX3fuel = HX3duty*1.3; %MW
HX4opcost = HX4duty*Price.Fuel*0.0075*5; %$
HX5opcost = HX5duty*Price.Fuel*0.0075; %$
HX6opcost = HX6duty*Price.Fuel*0.0075; %$

HXcost = sum([HX1cost,HX2cost,HX3cost,HX4cost,HX5cost,HX6cost]);
HXfuelreq = HX2fuel + HX3fuel; %MW
HXopcost = sum([HX4opcost,HX5opcost,HX6opcost]);


%reactor
RD = 1.09*3.28; %ft
RL = 10*3.28; %ft
RFc = 1*3.67;
Rcost = MS/280*101.9*RD^(1.066)*RL^(0.82)*(2.18+RFc); %$
Rduty = 0.585;
Rfuelreq = Rduty*1.13; %MW
RV = (((1.009/2)^2)*pi*10); %m^3
CatMass = RV*6350/1000; %MT

%sep system
%exchangers
D1condHXa = 335*10.7639; %ft^2
D1rebHXa = 260*10.7639; %ft^2
D2condHXa = 51*10.7639; %ft^2
D2rebHXa = 192*10.7639; %ft^2

D1condHXFc = (0.85+0)*3.75;
D1rebHXFc = (0.85+0)*3.75;
D2condHXFc = (0.85+0)*3.75;
D2rebHXFc = (0.85+0)*3.75;

D1condHXcost = MS/280*101.3*D1condHXa^(0.65)*(2.29+D1condHXFc); %$
D1rebHXcost = MS/280*101.3*D1rebHXa^(0.65)*(2.29+D1rebHXFc); %$
D2condHXcost = MS/280*101.3*D2condHXa^(0.65)*(2.29+D2condHXFc); %$
D2rebHXcost = MS/280*101.3*D2rebHXa^(0.65)*(2.29+D2rebHXFc); %$

%operation costs for sep HXers
D1rebHXduty = 9.167; %MW
D2rebHXduty = 5.131; %MW
D1condHXduty = 9.306*3600*8400/1055; %MM BTU/year
D2condHXduty = 4.222*3600*8400/1055; %MM BTU/year

D1rebHXfuel = D1rebHXduty*1.13; %MW
D2rebHXfuel = D2rebHXduty*1.13; %MW

Sepfuelreq = sum([D1rebHXfuel,D2rebHXfuel]); %MW

D1condHXopcost = D1condHXduty*Price.Fuel*0.0075; %$
D2condHXopcost = D2condHXduty*Price.Fuel*0.0075; %$

%columns
P = 1;
Fm = 3.5;
P0 = 4.5; %bar
t_p = 0.13*(P-P0)./P0;
Fp = 1+t_p*(1+exp(-t_p./2));
Fi = 1.38;
Fd = 3;
Fs = 1; % (for 60 cm tray spacing)
Ft = 0; %sieve
d0 = 1; %m
H0 = 6.1; %m
as = 0.82;
at = 1.8;
Fcshell = (Fm.*Fp -1+Fi.*Fd);
Fctray = (Fs + Ft + Fm);

D1d = mean([2.405,2.625]); %m (avg of stripping and rectifying)
D1H = 13.41+14.63; %m
D1V = (pi*(D1d/2)^2)*D1H;
CostShell = (MS./301).*(Fcshell).*5000.*(D1d./d0)*((D1H./H0).^as); %$
CostTray = (MS./301).*(Fctray).*5000.*((D1d./d0).^at).*(D1H./H0); %$
D1ColCost = CostShell + CostTray; %$

P = 20;
Fm = 3.5;
P0 = 4.5; %bar
t_p = 0.13*(P-P0)./P0;
Fp = 1+t_p*(1+exp(-t_p./2));
Fi = 1.38;
Fd = 3;
Fs = 1; % (for 60 cm tray spacing)
Ft = 0; %sieve
d0 = 1; %m
H0 = 6.1; %m
as = 0.82;
at = 1.8;
Fcshell = (Fm.*Fp -1+Fi.*Fd);
Fctray = (Fs + Ft + Fm);

D2d = mean([3.464,3.561]); %m (avg of stripping and rectifying)
D2H = 17.07+16.46; %m
D2V = (pi*(D2d/2)^2)*D2H;
CostShell = (MS./301).*(Fcshell).*5000.*(D2d./d0)*((D2H./H0).^as); %$
CostTray = (MS./301).*(Fctray).*5000.*((D2d./d0).^at).*(D2H./H0); %$
D2ColCost = CostShell + CostTray; %$


Sepcost = sum([D1condHXcost,D1rebHXcost,D2condHXcost,D2rebHXcost,...
    D1ColCost,D2ColCost]);
SepOpcost = sum([D1condHXopcost,D2condHXopcost]);

Fuel = (HXfuelreq+Sepfuelreq+Rfuelreq)*3600*8400/1055; %MM BTU/yr
FuelRefund = (DEflow*2726.3*74.123/1e6 + Hflow*286*2.016/1e6); %kJ/year
FuelRefund = 0.7*FuelRefund/1055000; %MM BTU/yr
% Fuel = Fuel-FuelRefund;
CO2prod = 4*(DEflow/74.123)*44.01; %MT/yr
CO2perMassProd = CO2prod/EAflow; %kg/kg

%fixed Costs
ISBL = Sepcost + Rcost + HXcost + 160000 + 100000; %$
FCI = 2.28*ISBL; %$
FCI = FCI*FCImultiplier;
SU = 0.1*FCI; %$
SU = SU*SUmultiplier;

Salvage = 0.01*FCI; %$
WC = Ethflow*Price.Eth*1/12; %$/yr
WC = WC*WCmultiplier;


TI = FCI + SU + WC; %$
FuelCost = Price.Fuel*Fuel; %$/yr
OPEX = FuelCost + SepOpcost + HXopcost + 70000 + 20000; %$
Revenue = EAflow*Price.EA; %$/yr
AGS = 0.1*Revenue; %$
ReactantCost = Ethflow*Price.Eth; %$/yr
CatalystCost = CatMass*Price.Cat; %$/yr
ConstructionLoan = FCI + WC + SU; %$
CO2cost = CO2prod*Price.CO2; %$/yr
COM = OPEX + AGS + CO2cost + ReactantCost; %$
FCneg3 = -0*FCI;
FCneg2 = -0*FCI;
FCneg1 = -0.6*FCI;
FC0 = -0.4*FCI;
WC0 = -WC;
SU0 = -SU;

FCneg3disc = (1/(1+InterestRate)^(-3))*FCneg3;
FCneg2disc = (1/(1+InterestRate)^(-2))*FCneg2;
FCneg1disc = (1/(1+InterestRate)^(-1))*FCneg1;
FC0disc = FC0;
WC0disc = WC0;
SU0disc = SU0;

BondPrinciple = sum([FCneg3disc,FCneg2disc,FCneg1disc,FC0disc,WC0disc,SU0disc]);
BondInterestPayments = BondPrinciple*(InterestRate-0.02);

ProfitBeforeTax = Revenue - COM;
Depreciation = (FCI + SU - Salvage)/10; %$/yr

for i = 1:15
    if i == 1
        Productivity = 0.5;
    else
        Productivity = 1;
    end

        CatYear = 1;
    
    if i<=10
        DepYear = 1;
    else
        DepYear = 0;
    end
    BondFinancing(i) = BondInterestPayments*DepYear;
    ProfitBT(i) = Productivity*ProfitBeforeTax - CatalystCost*CatYear*2 + BondFinancing(i);
    DepreciationAllowed(i) = Depreciation*DepYear;
    ProfitAfterTaxes(i) = (ProfitBT(i) - ...
        DepreciationAllowed(i))*(1-IncomeTaxRate);
    cashflow(i) = ProfitAfterTaxes(i) + DepreciationAllowed(i);
    
end

cashflow(10) = cashflow(10) + BondPrinciple + ((ProfitBT(i) - ...
    DepreciationAllowed(i))*(IncomeTaxRate));
cashflow(15) = cashflow(15) + WC + Salvage;
for i = 1:15
    disccashflow(i) = (1/(1+EnterpriseRate)^i)*cashflow(i);
end

NPV0 = sum(disccashflow);
NPVproject = NPV0*(1+EnterpriseRate)^(-2); %$
NPV = NPVproject/1e6;

end
