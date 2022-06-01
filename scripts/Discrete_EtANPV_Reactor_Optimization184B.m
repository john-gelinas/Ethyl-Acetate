clear
close all
tic
%this script solves the concentration profiles in our isothermal ideal pfr
%neglecting pressure drop to determine the selectivities, toluene
%conversion, and molar flow streams for a desired production rate of
%para-xylene by toluene disproportionation

%additionally this script determines revenue from flow streams from the
%separated components, and calculates NPV, TCI, ROI, and IRR from cash
%flows
Price.EtA = 1100;%[$/MT]
Price.Ethanol = 500;%[$/MT]
Price.Fuel = 2.5;%[$/MM BTU]
Price.Electricity = 0.08; % $/kw hr
Price.CO2 = 40; %$/MT
Price.Cat = 10000; %$/MT
MW.EtA = 88.11;%[g/gmol]
MW.Ethanol = 46.07;%[g/gmol]
MW.Acetaldehyde = 44.05;%[g/gmol]
MW.Ether = 74.12;%[g/gmol]
MW.H2 = 2.02;%[g/gmol]

lambda.EtA = 38.0; %kJ/mol, latent heat of evaporation
lambda.Ethanol = 42.3; %kJ/mol, latent heat of evaporation
lambda.Ether = 33.8; %kJ/mol, latent heat of evaporation

Cp.EthanolGas = 74/1000; %kJ/mol*K
Cp.EthanolLiq = 118/1000; %kJ/mol*K

dHcomb.H2 = 286; % kJ/mol

MethaneEnergy = 891; %kJ/gmol
MethaneCarbon = 1; %kJ/gmol

EnterpriseRate = 0.1;
ConstructionInterestRate = 0.05;
IncomeTaxRate = 0.27;
FCImultiplier = 1;
SUmultiplier = 1;
WCmultiplier = 1;

% P = 10; %[atm]
R = 0.08206; %[L*atm/K*mol]
R = R*1.01325; %[L*bar/K*mol]
F = 50;%[gmol/hr] basis flowrate used to solve system before scale up
PEtA = 67556; %[gmol/hr]
Vmax = 10; %[L]
labelT = {'T225' 'T230' 'T235' 'T240'};
labelP = {'P2' 'P4' 'P6' 'P8' 'P10'};
cc=lines(length(labelP));
for j = 2:2:10
    P = j; %bar
    countP = (j-2)/2 + 1;
    for i = 225:5:240
        countT = (i - 225)/5 + 1;
        T = i + 273;
        %assume reactor feed is pure toluene and equal to plant feed stream of
        %toluene
        Q = [1 0 0 0 0 0]*F; %feed molar flowrate in gmol/hr
        %solve ode
        resolution = 500;
        Vspan = [0 1.87 2.68 3.85 5.7 9.4];
        %         Vspan = logspace(-3,log10(Vmax),1001);
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [Vout, Nout] = ode15s(@(V,N) EtAcetateReactor_PBR_Ideal_184B(V,N,T,P),Vspan,Q,opts);
        %store output of ode45 in structures containing all values at different
        %temperatures
        V.(labelP{countP}).(labelT{countT}) = Vout(2:end);
        Ethanol.(labelP{countP}).(labelT{countT}) = Nout(2:end,1);
        Acetaldyhyde.(labelP{countP}).(labelT{countT}) = Nout(2:end,2);
        H2.(labelP{countP}).(labelT{countT}) = Nout(2:end,3);
        EtA.(labelP{countP}).(labelT{countT}) = Nout(2:end,4);
        Ether.(labelP{countP}).(labelT{countT}) = Nout(2:end,5);
        H2O.(labelP{countP}).(labelT{countT}) = Nout(2:end,6);
        
        %compute tau for all reactor sizes output from ode45
        for new = 1:length(V.(labelP{countP}).(labelT{countT}))
            tau.(labelP{countP}).(labelT{countT})(new,1) = Vout(new)*P./(R*T*sum(Nout(new,:)));
        end
        %calculate selectivity and conversion for selected T and Ft
        Ft = F;
        xt.(labelP{countP}).(labelT{countT}) = (Ft - Ethanol.(labelP{countP}).(labelT{countT}))./Ft;%reactor conversion of Ethanol
        sEtA.(labelP{countP}).(labelT{countT})  = EtA.(labelP{countP}).(labelT{countT})./((1/1)*(Ft - Ethanol.(labelP{countP}).(labelT{countT})));%reactor selectivity for EtA
        sEther.(labelP{countP}).(labelT{countT}) = Ether.(labelP{countP}).(labelT{countT})./((1/1)*(Ft - Ethanol.(labelP{countP}).(labelT{countT})));%reactor selectivity for Ether
        
        %calculate real molar feeds, volumetric flow, and total reactor volume
        %at full scale
        FEthanol.(labelP{countP}).(labelT{countT}) = PEtA./(xt.(labelP{countP}).(labelT{countT}).*sEtA.(labelP{countP}).(labelT{countT}));%feed reate of Ethanol to reactor for desired PEtA [gmol/hr]
        qreal.(labelP{countP}).(labelT{countT}) = R*T*FEthanol.(labelP{countP}).(labelT{countT})/P;%real volumetric flowrate [L/hr]
        Vreal.(labelP{countP}).(labelT{countT}) = tau.(labelP{countP}).(labelT{countT}).*qreal.(labelP{countP}).(labelT{countT});%real reactor total volume [L]
        
        %calculate the product molar flowrates from full scale reactor
        PEther.(labelP{countP}).(labelT{countT}) = FEthanol.(labelP{countP}).(labelT{countT}).*xt.(labelP{countP}).(labelT{countT}).*sEther.(labelP{countP}).(labelT{countT});%Product molar flow Ether [gmol/hr]
        PH2O.(labelP{countP}).(labelT{countT}) = PEther.(labelP{countP}).(labelT{countT});
        PEthanol.(labelP{countP}).(labelT{countT}) = FEthanol.(labelP{countP}).(labelT{countT}).*(1-xt.(labelP{countP}).(labelT{countT}));%Product molar flow Ethanol [gmol/hr]
        PH2.(labelP{countP}).(labelT{countT}) = FEthanol.(labelP{countP}).(labelT{countT}) - PEthanol.(labelP{countP}).(labelT{countT});%Product molar flow Hydrogen [gmol/hr]
        PTotal.(labelP{countP}).(labelT{countT}) = PEther.(labelP{countP}).(labelT{countT})+PEthanol.(labelP{countP}).(labelT{countT})+PH2.(labelP{countP}).(labelT{countT})+PH2O.(labelP{countP}).(labelT{countT})+PEtA;
        
        %create vectors of mass flowrates
        PMEther.(labelP{countP}).(labelT{countT}) = PEther.(labelP{countP}).(labelT{countT}).*MW.Ether./(10^3);%[kg/hr]
        PMEthanol.(labelP{countP}).(labelT{countT}) = PEthanol.(labelP{countP}).(labelT{countT}).*MW.Ethanol./(10^3);%[kg/hr]
        FMEthanol.(labelP{countP}).(labelT{countT}) = FEthanol.(labelP{countP}).(labelT{countT}).*MW.Ethanol./(10^3);%[kg/hr]
        PMEtA.(labelP{countP}).(labelT{countT}) = (PEtA*MW.EtA/(10^3)).*ones(length(PEthanol.(labelP{countP}).(labelT{countT})),1);%[kg/hr]
        
        %calculate the fresh feed stream flowrate assuming perfect separations
        FreshFeed.(labelP{countP}).(labelT{countT}) = FEthanol.(labelP{countP}).(labelT{countT})-PEthanol.(labelP{countP}).(labelT{countT});
        FreshFeedM.(labelP{countP}).(labelT{countT}) = FreshFeed.(labelP{countP}).(labelT{countT}).*MW.Ethanol./(10^3);%[kg/hr]
        
        %%equipment costs
        %reactor - sufficiently estimated with single tubular reactor cost
        VR.(labelP{countP}).(labelT{countT}) = Vreal.(labelP{countP}).(labelT{countT})/1000; %m^3
        VR.(labelP{countP}).(labelT{countT}) = VR.(labelP{countP}).(labelT{countT})*(3.28^3); %ft^3
        Fm = 2.25; %stainless steel - H2 embrittlement
        % pressure correction
        if P < 3.45
            Fp = 1.00;
        elseif P < 6.9
            Fp = 1.05;
        elseif P < 13.8
            Fp = 1.15;
        elseif P < 20.7
            Fp = 1.20;
        elseif P < 27.6
            Fp = 1.35;
        elseif P < 34.5
            Fp = 1.45;
        end
        Fc = Fm*Fp;
        MS = 1650;
        %set L:D ratio
        D.(labelP{countP}).(labelT{countT}) = ((4/10)*(VR.(labelP{countP}).(labelT{countT})/pi)).^(1/3); %ft
        H.(labelP{countP}).(labelT{countT}) = 10*D.(labelP{countP}).(labelT{countT}); %ft
        Area.(labelP{countP}).(labelT{countT}) = pi.*(D.(labelP{countP}).(labelT{countT})./2).^2; %ft^2
        ReactorCost.(labelP{countP}).(labelT{countT}) = (MS/280)*101.9*(D.(labelP{countP}).(labelT{countT}).^1.066).*(H.(labelP{countP}).(labelT{countT}).^0.82)*(2.18+Fc);
        ReactorCost.(labelP{countP}).(labelT{countT}) = ReactorCost.(labelP{countP}).(labelT{countT})/1e6; %MM$/yr
        
        
        % feed furnace
        eff = 0.70; % typical furnace efficiency
        bp.Ethanol = 8*P+70; %linear approximation from phase diagram
        Tf1 = 25;
        Tf2 = bp.Ethanol;
        T01 = bp.Ethanol;
        T02 = T - 273;
        QH2.(labelP{countP}).(labelT{countT}) = (eff)*PH2.(labelP{countP}).(labelT{countT})*dHcomb.H2; %kJ/hr
        Qtotal.(labelP{countP}).(labelT{countT}) = FEthanol.(labelP{countP}).(labelT{countT})*(lambda.Ethanol + ...
            (Tf2-Tf1)*Cp.EthanolLiq + (T02-T01)*Cp.EthanolGas); %kJ/hr
        for c = 1:length(QH2.(labelP{countP}).(labelT{countT}))
            if QH2.(labelP{countP}).(labelT{countT})(c) > Qtotal.(labelP{countP}).(labelT{countT})(c)
                Qfurnace.(labelP{countP}).(labelT{countT})(c,1) = Qtotal.(labelP{countP}).(labelT{countT})(c);
                Qresid.(labelP{countP}).(labelT{countT})(c,1) = 0;
            elseif Qtotal.(labelP{countP}).(labelT{countT})(c) > QH2.(labelP{countP}).(labelT{countT})(c)
                Qfurnace.(labelP{countP}).(labelT{countT})(c,1) = QH2.(labelP{countP}).(labelT{countT})(c);
                Qresid.(labelP{countP}).(labelT{countT})(c,1) = Qtotal.(labelP{countP}).(labelT{countT})(c) - QH2.(labelP{countP}).(labelT{countT})(c);
            end
        end
        Qfurnace.(labelP{countP}).(labelT{countT}) = Qfurnace.(labelP{countP}).(labelT{countT})*0.947817/1e6; %MM BTU/hr
        %         FurnaceCarbon.(labelP{countP}).(labelT{countT}) = H2carbon*PH2.(labelP{countP}).(labelT{countT}); %gmol/hr
        %         FurnaceCarbon.(labelP{countP}).(labelT{countT}) = FurnaceCarbon.(labelP{countP}).(labelT{countT})*44.01/1000; %kg/hr
        %         FurnaceCarbon.(labelP{countP}).(labelT{countT}) = FurnaceCarbon.(labelP{countP}).(labelT{countT})*(24*50*7)/1000; %MT/year
        FurnaceCarbon.(labelP{countP}).(labelT{countT}) = zeros(length(QH2.(labelP{countP}).(labelT{countT})),1); %no carbon production from H2 combustion
        Fc = 1.5; %stainless steel
        FurnaceCost.(labelP{countP}).(labelT{countT}) = (MS/280)*(5.07e3)...
            *(Qfurnace.(labelP{countP}).(labelT{countT}).^0.85)*(1.23+Fc); %$
        FurnaceCost.(labelP{countP}).(labelT{countT}) = FurnaceCost.(labelP{countP}).(labelT{countT})/1e6; %MM$
        
        
        % feed heat exchangers
        SteamTemp = 139; %50psig sat steam
        
        for c = 1:length(QH2.(labelP{countP}).(labelT{countT}))
            if QH2.(labelP{countP}).(labelT{countT})(c) > Qtotal.(labelP{countP}).(labelT{countT})(c)
                %heat exchanger 1:
                T01 = 25; %C
                Tf1(c,1) = Qresid.(labelP{countP}).(labelT{countT})(c)./(FEthanol.(labelP{countP}).(labelT{countT})(c).*Cp.EthanolLiq);
                dTA = SteamTemp - T01; %inlet temp difference
                dTB(c,1) = SteamTemp - Tf1(c); %outlet temp difference
                %sensible heat
                FeedHeating.(labelP{countP}).(labelT{countT})(c,1) = Qresid.(labelP{countP}).(labelT{countT})(c); %kJ/hr
                FeedHeatinghr.(labelP{countP}).(labelT{countT})(c,1) = FeedHeating.(labelP{countP}).(labelT{countT})(c)*0.947817/1e6; %MM Btu/hr
                H1SteamCost.(labelP{countP}).(labelT{countT})(c,1) = FeedHeatinghr.(labelP{countP}).(labelT{countT})(c)*Price.Fuel*0.7; %$/hr (from table E.1-1)
                H1SteamCost.(labelP{countP}).(labelT{countT})(c,1) = H1SteamCost.(labelP{countP}).(labelT{countT})(c)*(24*7*50)/1e6; %MM$/yr
                H1SteamCost.(labelP{countP}).(labelT{countT})(c,1) = 0; %MM$/yr % assume heat integration with feed effluent
                % will be ingnored as design includes HX integration
                dTlm(c,1) = (dTA - dTB(c))./(log(dTA./dTB(c))); %approximation for log mean temp diff
                dTlm(c,1) = dTlm(c)*9/5 + 32; %F
            elseif Qtotal.(labelP{countP}).(labelT{countT})(c) > QH2.(labelP{countP}).(labelT{countT})(c)
                dTlm(c,1) = 1;
                FeedHeatinghr.(labelP{countP}).(labelT{countT})(c,1) = 0;
                H1SteamCost.(labelP{countP}).(labelT{countT})(c,1) = 0; %MM$/yr % assume heat integration with feed effluent
            end
        end
        U = 150; %BTU/ft^2 F h
        eta = 0.7;
        A = (eta)*FeedHeatinghr.(labelP{countP}).(labelT{countT})*1e6./(U.*dTlm); %ft^2
        Fm = 2.81; %material correction, stainless steel
        Fd = 1; %design correction
        if P < 20.7
            Fp = 0; %pressure correction
        elseif P < 27.6
            Fp = 0.1;
        elseif P < 55.2
            Fp = 0.25;
        end
        Fc = (Fd + Fp)*Fm;
        HX1Q.(labelP{countP}).(labelT{countT}) = FeedHeating.(labelP{countP}).(labelT{countT}); %kJ/hr
        HXcost1.(labelP{countP}).(labelT{countT}) = (MS/280)*101.3*(A.^0.65)*(2.29+Fc); %$
        
        HXcost.(labelP{countP}).(labelT{countT}) = HXcost1.(labelP{countP}).(labelT{countT})/1e6; %MM$
        
        HeatOPCost.(labelP{countP}).(labelT{countT}) = H1SteamCost.(labelP{countP}).(labelT{countT}); %MM$ /yr
        
        %%separation system
        PSepTotal.(labelP{countP}).(labelT{countT}) = PEthanol.(labelP{countP}).(labelT{countT})+PH2O.(labelP{countP}).(labelT{countT})+PEtA;
        XH2O.(labelP{countP}).(labelT{countT}) = PH2O.(labelP{countP}).(labelT{countT})./PSepTotal.(labelP{countP}).(labelT{countT});
        XEthanol.(labelP{countP}).(labelT{countT}) = PEthanol.(labelP{countP}).(labelT{countT})./PSepTotal.(labelP{countP}).(labelT{countT});
        XEtA.(labelP{countP}).(labelT{countT}) = PEtA./PSepTotal.(labelP{countP}).(labelT{countT});
        
        for y = 1:length(xt.(labelP{countP}).(labelT{countT}))
            %Aspen Plus data
            if xt.(labelP{countP}).(labelT{countT})(y) > 0.47
                Dist1Qc.(labelP{countP}).(labelT{countT})(y) = 5.46*3.6e6; %kj/hr
                Dist1Qr.(labelP{countP}).(labelT{countT})(y) = 6.22*3.6e6; %kj/hr
                Dist1Nreal.(labelP{countP}).(labelT{countT})(y) = 38;
                Dist1rhog.(labelP{countP}).(labelT{countT})(y) = 2.36;
                Dist1V.(labelP{countP}).(labelT{countT})(y) = 572549;
                Dist1MV.(labelP{countP}).(labelT{countT})(y) = 57.01357;
                Dist2Qc.(labelP{countP}).(labelT{countT})(y) = 3.72*3.6e6; %kj/hr
                Dist2Qr.(labelP{countP}).(labelT{countT})(y) = 3.33*3.6e6; %kj/hr
                Dist2Nreal.(labelP{countP}).(labelT{countT})(y) = 48;
                Dist2rhog.(labelP{countP}).(labelT{countT})(y) = 2.41;
                Dist2V.(labelP{countP}).(labelT{countT})(y) = 561793;
                Dist2MV.(labelP{countP}).(labelT{countT})(y) = 66.80337;
            elseif xt.(labelP{countP}).(labelT{countT})(y) > 0.42
                Dist1Qc.(labelP{countP}).(labelT{countT})(y) = 8.4*3.6e6; %kj/hr
                Dist1Qr.(labelP{countP}).(labelT{countT})(y) = 9.58*3.6e6; %kj/hr
                Dist1Nreal.(labelP{countP}).(labelT{countT})(y) = 28;
                Dist1rhog.(labelP{countP}).(labelT{countT})(y) = 2.34;
                Dist1V.(labelP{countP}).(labelT{countT})(y) = 881195;
                Dist1MV.(labelP{countP}).(labelT{countT})(y) = 56.54484;
                Dist2Qc.(labelP{countP}).(labelT{countT})(y) = 3.66*3.6e6; %kj/hr
                Dist2Qr.(labelP{countP}).(labelT{countT})(y) = 3.28*3.6e6; %kj/hr
                Dist2Nreal.(labelP{countP}).(labelT{countT})(y) = 43;
                Dist2rhog.(labelP{countP}).(labelT{countT})(y) = 2.41;
                Dist2V.(labelP{countP}).(labelT{countT})(y) = 553393;
                Dist2MV.(labelP{countP}).(labelT{countT})(y) = 66.90638;
            elseif xt.(labelP{countP}).(labelT{countT})(y) > 0.37
                Dist1Qc.(labelP{countP}).(labelT{countT})(y) = 8.32*3.6e6; %kj/hr
                Dist1Qr.(labelP{countP}).(labelT{countT})(y) = 9.5*3.6e6; %kj/hr
                Dist1Nreal.(labelP{countP}).(labelT{countT})(y) = 33;
                Dist1rhog.(labelP{countP}).(labelT{countT})(y) = 2.33;
                Dist1V.(labelP{countP}).(labelT{countT})(y) = 876150;
                Dist1MV.(labelP{countP}).(labelT{countT})(y) = 56.14021;
                Dist2Qc.(labelP{countP}).(labelT{countT})(y) = 3.6*3.6e6; %kj/hr
                Dist2Qr.(labelP{countP}).(labelT{countT})(y) = 3.25*3.6e6; %kj/hr
                Dist2Nreal.(labelP{countP}).(labelT{countT})(y) = 42;
                Dist2rhog.(labelP{countP}).(labelT{countT})(y) = 2.4;
                Dist2V.(labelP{countP}).(labelT{countT})(y) = 547246;
                Dist2MV.(labelP{countP}).(labelT{countT})(y) = 66.729;
            elseif xt.(labelP{countP}).(labelT{countT})(y) > 0.32
                Dist1Qc.(labelP{countP}).(labelT{countT})(y) = 8.28*3.6e6; %kj/hr
                Dist1Qr.(labelP{countP}).(labelT{countT})(y) = 9.45*3.6e6; %kj/hr
                Dist1Nreal.(labelP{countP}).(labelT{countT})(y) = 47;
                Dist1rhog.(labelP{countP}).(labelT{countT})(y) = 2.31;
                Dist1V.(labelP{countP}).(labelT{countT})(y) = 869659;
                Dist1MV.(labelP{countP}).(labelT{countT})(y) = 56.24488;
                Dist2Qc.(labelP{countP}).(labelT{countT})(y) = 3.99*3.6e6; %kj/hr
                Dist2Qr.(labelP{countP}).(labelT{countT})(y) = 3.52*3.6e6; %kj/hr
                Dist2Nreal.(labelP{countP}).(labelT{countT})(y) = 37;
                Dist2rhog.(labelP{countP}).(labelT{countT})(y) = 2.4;
                Dist2V.(labelP{countP}).(labelT{countT})(y) = 598337.6;
                Dist2MV.(labelP{countP}).(labelT{countT})(y) = 67.65983;
            elseif xt.(labelP{countP}).(labelT{countT})(y) > 0
                Dist1Qc.(labelP{countP}).(labelT{countT})(y) = 18.6*3.6e6; %kj/hr
                Dist1Qr.(labelP{countP}).(labelT{countT})(y) = 21.2*3.6e6; %kj/hr
                Dist1Nreal.(labelP{countP}).(labelT{countT})(y) = 22;
                Dist1rhog.(labelP{countP}).(labelT{countT})(y) = 2.29;
                Dist1V.(labelP{countP}).(labelT{countT})(y) = 1947300;
                Dist1MV.(labelP{countP}).(labelT{countT})(y) = 55.09084;
                Dist2Qc.(labelP{countP}).(labelT{countT})(y) = 5.51*3.6e6; %kj/hr
                Dist2Qr.(labelP{countP}).(labelT{countT})(y) = 4.57*3.6e6; %kj/hr
                Dist2Nreal.(labelP{countP}).(labelT{countT})(y) = 55;
                Dist2rhog.(labelP{countP}).(labelT{countT})(y) = 2.39;
                Dist2V.(labelP{countP}).(labelT{countT})(y) = 796140;
                Dist2MV.(labelP{countP}).(labelT{countT})(y) = 66.729;

            end
            
            %distillation column 1
            Dist1P = 1; %bar

%             Dist1MvFeed = xa.*MW.Ether + xb.*MW.Ethanol + xc.*MW.EtA + xd.*MW.Ether;
%             Dist1MvDist = xa.*MW.Ether./xa;
%             Dist1MvBottoms = (xb.*MW.Ethanol + xc.*MW.EtA + xd.*MW.Ether).*(1./(xb+xc+xd));
%             Dist1lambdaD = lambda.Benzene; %latent heat of distillate mixture
%             Dist1lambdaB = (xb.*lambda.Tol + xc.*lambda.ParaX + xd.*lambda.TMB).*(1./(xb+xc+xd)); %latent heat of bottoms mixture

            Dist1DSattemp = 85; %deg C, BP of mixture estimate
            Dist1BSattemp = 76.5; %deg C, BP of Ethyl acetate
            Dist1CoolingWaterCoolTemp = 25; % deg C
            Dist1CoolingWaterHotTemp = Dist1CoolingWaterCoolTemp + 15; % deg C
            Dist1SteamTemp = 194; % deg C, 200 psia steam
            Dist1dTlmD.(labelP{countP}).(labelT{countT}) = ...
                ((Dist1DSattemp - Dist1CoolingWaterCoolTemp) - (Dist1DSattemp - Dist1CoolingWaterHotTemp))/...
                log((Dist1DSattemp - Dist1CoolingWaterCoolTemp)/(Dist1DSattemp - Dist1CoolingWaterHotTemp));
            Dist1dTlmB.(labelP{countP}).(labelT{countT}) = (Dist1SteamTemp - Dist1BSattemp);
            Dist1UD = 800; %W/m^2 K
            Dist1UD = Dist1UD*60*60/1000; %kJ/hr m^2 K
            Dist1UB = 820; %W/m^2 K
            Dist1UB = Dist1UB*60*60/1000; %kJ/hr m^2 K
            Dist1AreaD.(labelP{countP}).(labelT{countT}) = Dist1Qc.(labelP{countP}).(labelT{countT})./(Dist1UD.*Dist1dTlmD.(labelP{countP}).(labelT{countT})); %m^2
            Dist1AreaB.(labelP{countP}).(labelT{countT}) = Dist1Qr.(labelP{countP}).(labelT{countT})./(Dist1UB.*Dist1dTlmB.(labelP{countP}).(labelT{countT})); %m^2
            
            Fd = 1;
            P0 = 11; %bar
            Fp = 0.10*(Dist1P-P0)./P0;
            Fm = 1;
            Fi = 1.38;
            FD = 2.3;
            ah = 0.65;
            A0 = 93; %m^2
            Dist1HXDCost.(labelP{countP}).(labelT{countT}) = ... %$
                (MS./301).*((Fd+Fp).*Fm-1+Fi.*FD).*8700.*((Dist1AreaD.(labelP{countP}).(labelT{countT})./A0).^ah);
            Dist1HXBCost.(labelP{countP}).(labelT{countT}) = ... %$
                (MS./301).*((Fd+Fp).*Fm-1+Fi.*FD).*8700.*((Dist1AreaB.(labelP{countP}).(labelT{countT})./A0).^ah);
            
            %column
            Ht = 60; %cm
            Hmin = 3.*Ht; %cm
            Dist1H.(labelP{countP}).(labelT{countT}) = Hmin + Ht.*Dist1Nreal.(labelP{countP}).(labelT{countT}); %cm
            Dist1H.(labelP{countP}).(labelT{countT}) = Dist1H.(labelP{countP}).(labelT{countT})./100; %m
            c0 = 439;
            phi = 0.6;
            rhol = 870; %kg/m^3
            rhog = Dist1rhog.(labelP{countP}).(labelT{countT}); %kg/m^3
            DV = Dist1V.(labelP{countP}).(labelT{countT});
            Mv = Dist1MV.(labelP{countP}).(labelT{countT});
            DistArea = (Mv./sqrt(rhol.*rhog))*(1./(phi.*c0)).*(.1./0.8).*DV;
            Diameter = 2.*sqrt(DistArea./pi);
            Fm = 1;
            P0 = 4.5; %bar
            t_p = 0.13*(Dist1P-P0)./P0;
            Fp = 1+t_p*(1+exp(-t_p./2));
            Fi = 1.38;
            Fd = 3;
            Fs = 1; % (for 60 cm tray spacing)
            Ft = 0; %sieve
            Dist1D.(labelP{countP}).(labelT{countT}) = Diameter;
            d0 = 1; %m
            H0 = 6.1; %m
            as = 0.82;
            at = 1.8;
            CostShell = (MS./301).*(Fm.*Fp -1+Fi.*Fd).*5000.*(Dist1D.(labelP{countP}).(labelT{countT})./d0).*((Dist1H.(labelP{countP}).(labelT{countT})./H0).^as); %$
            CostTray = (MS./301).*(Fs + Ft + Fm).*5000.*((Dist1D.(labelP{countP}).(labelT{countT})./d0).^at).*(Dist1H.(labelP{countP}).(labelT{countT})./H0); %$
            Dist1ColumnCost.(labelP{countP}).(labelT{countT}) = (CostShell + CostTray)/10; %$
            
            Dist1QcBTU.(labelP{countP}).(labelT{countT}) = Dist1Qc.(labelP{countP}).(labelT{countT})*0.947817/1e6; %MM Btu/hr
            Dist1HXDOpCost.(labelP{countP}).(labelT{countT}) = Dist1QcBTU.(labelP{countP}).(labelT{countT})*Price.Fuel*1.13; %$/hr
            Dist1HXDOpCost.(labelP{countP}).(labelT{countT}) = Dist1HXDOpCost.(labelP{countP}).(labelT{countT})*(24*7*50)/1e6; %MM$/yr
            
            Dist1QrBTU.(labelP{countP}).(labelT{countT}) = Dist1Qr.(labelP{countP}).(labelT{countT})*0.947817/1e6; %MM Btu/hr
            Dist1HXBOpCost.(labelP{countP}).(labelT{countT}) = Dist1QrBTU.(labelP{countP}).(labelT{countT})*Price.Fuel*0.93; %$/hr
            Dist1HXBOpCost.(labelP{countP}).(labelT{countT}) = Dist1HXBOpCost.(labelP{countP}).(labelT{countT})*(24*7*50)/1e6; %MM$/yr
            
            %distillation column 2
            Dist2P = 20; %bar

%             Dist2MvFeed = xa.*MW.Ethanol + xb.*MW.EtA + xc.*MW.Ether;
%             Dist2MvDist = xa.*MW.Ethanol./xa;
%             Dist2MvBottoms = (xb.*MW.EtA + xc.*MW.Ether).*(1./(xb+xc));

%             Dist2lambdaD = lambda.Tol; %latent heat of distillate mixture
%             Dist2lambdaB = (xb.*lambda.ParaX + xc.*lambda.TMB).*(1./(xb+xc)); %latent heat of bottoms mixture

            Dist2DSattemp = 80; %deg C, BP of mixture, estimate
            Dist2BSattemp = 90; %deg C, BP of bottoms mixture, estimate
            Dist2CoolingWaterCoolTemp = 25; % deg C
            Dist2CoolingWaterHotTemp = Dist1CoolingWaterCoolTemp + 15; % deg C
            Dist2SteamTemp = 194; % deg C, 200 psia steam
            Dist2dTlmD.(labelP{countP}).(labelT{countT}) = ...
                ((Dist2DSattemp - Dist2CoolingWaterCoolTemp) - (Dist2DSattemp - Dist2CoolingWaterHotTemp))/...
                log((Dist2DSattemp - Dist2CoolingWaterCoolTemp)/(Dist2DSattemp - Dist2CoolingWaterHotTemp));
            Dist2dTlmB.(labelP{countP}).(labelT{countT}) = (Dist2SteamTemp - Dist2BSattemp);
            Dist2UD = 800; %W/m^2 K
            Dist2UD = Dist2UD*60*60/1000; %kJ/hr m^2 K
            Dist2UB = 820; %W/m^2 K
            Dist2UB = Dist2UB*60*60/1000; %kJ/hr m^2 K
            Dist2AreaD.(labelP{countP}).(labelT{countT}) = Dist2Qc.(labelP{countP}).(labelT{countT})./(Dist2UD.*Dist2dTlmD.(labelP{countP}).(labelT{countT})); %m^2
            Dist2AreaB.(labelP{countP}).(labelT{countT}) = Dist2Qr.(labelP{countP}).(labelT{countT})./(Dist2UB.*Dist2dTlmB.(labelP{countP}).(labelT{countT})); %m^2
            
            Fd = 1;
            P0 = 11; %bar
            Fp = 0.10*(Dist2P-P0)./P0;
            Fm = 1;
            Fi = 1.38;
            FD = 2.3;
            ah = 0.65;
            A0 = 93; %m^2
            Dist2HXDCost.(labelP{countP}).(labelT{countT}) = ... %$
                (MS./301).*((Fd+Fp).*Fm-1+Fi.*FD).*8700.*((Dist2AreaD.(labelP{countP}).(labelT{countT})./A0).^ah);
            Dist2HXBCost.(labelP{countP}).(labelT{countT}) = ... %$
                (MS./301).*((Fd+Fp).*Fm-1+Fi.*FD).*8700.*((Dist2AreaB.(labelP{countP}).(labelT{countT})./A0).^ah);
            
            
            
            Ht = 60; %cm
            Hmin = 3.*Ht; %cm
            Dist2H.(labelP{countP}).(labelT{countT}) = Hmin + Ht.*Dist1Nreal.(labelP{countP}).(labelT{countT}); %cm
            Dist2H.(labelP{countP}).(labelT{countT}) = Dist2H.(labelP{countP}).(labelT{countT})./100; %m
            c0 = 439;
            phi = 0.6;
            rhol = 870; %kg/m^3
            rhog = Dist2rhog.(labelP{countP}).(labelT{countT}); %kg/m^3
            DV = Dist2V.(labelP{countP}).(labelT{countT});
            Mv = Dist2MV.(labelP{countP}).(labelT{countT});
            DistArea = (Mv./sqrt(rhol.*rhog))*(1./(phi.*c0)).*(.1./0.8).*DV;
            Diameter = 2.*sqrt(DistArea./pi);
            Fm = 1;
            P0 = 4.5; %bar
            t_p = 0.13*(Dist2P-P0)./P0;
            Fp = 1+t_p*(1+exp(-t_p./2));
            Fi = 1.38;
            Fd = 3;
            Fs = 1; % (for 60 cm tray spacing)
            Ft = 0; %sieve
            Dist2D.(labelP{countP}).(labelT{countT}) = Diameter;
            d0 = 1; %m
            H0 = 6.1; %m
            as = 0.82;
            at = 1.8;
            CostShell = (MS./301).*(Fm.*Fp -1+Fi.*Fd).*5000.*(Dist2D.(labelP{countP}).(labelT{countT})./d0).*((Dist2H.(labelP{countP}).(labelT{countT})./H0).^as); %$
            CostTray = (MS./301).*(Fs + Ft + Fm).*5000.*((Dist2D.(labelP{countP}).(labelT{countT})./d0).^at).*(Dist2H.(labelP{countP}).(labelT{countT})./H0); %$
            Dist2ColumnCost.(labelP{countP}).(labelT{countT}) = (CostShell + CostTray)/10; %$
            
            Dist2QcBTU.(labelP{countP}).(labelT{countT}) = Dist2Qc.(labelP{countP}).(labelT{countT})*0.947817/1e6; %MM Btu/hr
            Dist2HXDOpCost.(labelP{countP}).(labelT{countT}) = Dist2QcBTU.(labelP{countP}).(labelT{countT})*Price.Fuel*1.13; %$/hr
            Dist2HXDOpCost.(labelP{countP}).(labelT{countT}) = Dist2HXDOpCost.(labelP{countP}).(labelT{countT})*(24*7*50)/1e6; %MM$/yr
            
            Dist2QrBTU.(labelP{countP}).(labelT{countT}) = Dist2Qr.(labelP{countP}).(labelT{countT})*0.947817/1e6; %MM Btu/hr
            Dist2HXBOpCost.(labelP{countP}).(labelT{countT}) = Dist2QrBTU.(labelP{countP}).(labelT{countT})*Price.Fuel*0.93; %$/hr
            Dist2HXBOpCost.(labelP{countP}).(labelT{countT}) = Dist2HXBOpCost.(labelP{countP}).(labelT{countT})*(24*7*50)/1e6; %MM$/yr
            
%             %distillation column 3
%             Dist3P = 20; %bar
% 
% %             MvFeed = xa.*MW.EtA + xb.*MW.Ether;
% %             MvDist = xa.*MW.EtA./xa;
% %             MvBottoms = xb.*MW.Ether./xb;
% 
%             Dist3DSattemp = 78; %deg C, BP of Ethanol
%             Dist3BSattemp = 100; %deg C, BP of water
%             Dist3CoolingWaterCoolTemp = 25; % deg C
%             Dist3CoolingWaterHotTemp = Dist1CoolingWaterCoolTemp + 15; % deg C
%             Dist3SteamTemp = 194; % deg C, 200 psia steam
%             Dist1dTlmD.(labelP{countP}).(labelT{countT}) = ...
%                 ((Dist3DSattemp - Dist3CoolingWaterCoolTemp) - (Dist3DSattemp - Dist3CoolingWaterHotTemp))/...
%                 log((Dist3DSattemp - Dist3CoolingWaterCoolTemp)/(Dist3DSattemp - Dist3CoolingWaterHotTemp));
%             Dist3dTlmB.(labelP{countP}).(labelT{countT}) = (Dist3SteamTemp - Dist3BSattemp);
%             Dist3UD = 770; %W/m^2 K
%             Dist3UD = Dist3UD*60*60/1000; %kJ/hr m^2 K
%             Dist3UB = 820; %W/m^2 K
%             Dist3UB = Dist3UB*60*60/1000; %kJ/hr m^2 K
%             Dist3AreaD.(labelP{countP}).(labelT{countT}) = Dist3Qc.(labelP{countP}).(labelT{countT})./(Dist3UD.*Dist2dTlmD.(labelP{countP}).(labelT{countT})); %m^2
%             Dist3AreaB.(labelP{countP}).(labelT{countT}) = Dist3Qr.(labelP{countP}).(labelT{countT})./(Dist3UB.*Dist2dTlmB.(labelP{countP}).(labelT{countT})); %m^2
%             
%             Fd = 1;
%             P0 = 11; %bar
%             Fp = 0.10*(Dist3P-P0)./P0;
%             Fm = 1;
%             Fi = 1.38;
%             FD = 2.3;
%             ah = 0.65;
%             A0 = 93; %m^2
%             Dist3HXDCost.(labelP{countP}).(labelT{countT}) = ... %$
%                 (MS./301)*((Fd+Fp).*Fm-1+Fi.*FD).*8700.*((Dist3AreaD.(labelP{countP}).(labelT{countT})./A0).^ah);
%             Dist3HXBCost.(labelP{countP}).(labelT{countT}) = ... %$
%                 (MS./301)*((Fd+Fp).*Fm-1+Fi.*FD).*8700.*((Dist3AreaB.(labelP{countP}).(labelT{countT})./A0).^ah);
%             
%             
%             Ht = 60; %cm
%             Hmin = 3.*Ht; %cm
%             Dist3H.(labelP{countP}).(labelT{countT}) = Hmin + Ht.*Dist1Nreal.(labelP{countP}).(labelT{countT}); %cm
%             Dist3H.(labelP{countP}).(labelT{countT}) = Dist3H.(labelP{countP}).(labelT{countT})./100; %m
%             c0 = 439;
%             phi = 0.6;
%             rhol = 870; %kg/m^3
%             rhog = 1; %kg/m^3
%             Area = (Mv./sqrt(rhol.*rhog))*(1./(phi.*c0)).*(.1./0.8).*V;
%             Diameter = 2.*sqrt(Area./pi);
%             Fm = 1;
%             P0 = 4.5; %bar
%             t_p = 0.13*(Dist3P-P0)./P0;
%             Fp = 1+t_p*(1+exp(-t_p./2));
%             Fi = 1.38;
%             Fd = 3;
%             Fs = 1; % (for 60 cm tray spacing)
%             Ft = 0; %sieve
%             Dist3D.(labelP{countP}).(labelT{countT}) = Diameter;
%             d0 = 1; %m
%             H0 = 6.1; %m
%             as = 0.82;
%             at = 1.8;
%             CostShell = (MS./301).*(Fm.*Fp -1+Fi.*Fd).*5000.*(Dist3D.(labelP{countP}).(labelT{countT})./d0).*((Dist3H.(labelP{countP}).(labelT{countT})./H0).^as); %$
%             CostTray = (MS./301).*(Fs + Ft + Fm).*5000.*((Dist3D.(labelP{countP}).(labelT{countT})./d0).^at).*(Dist3H.(labelP{countP}).(labelT{countT})./H0); %$
%             Dist3ColumnCost.(labelP{countP}).(labelT{countT}) = (CostShell + CostTray)/10; %$
%             
%             Dist3QcBTU.(labelP{countP}).(labelT{countT}) = Dist3Qc.(labelP{countP}).(labelT{countT})*0.947817/1e6; %MM Btu/hr
%             Dist3HXDOpCost.(labelP{countP}).(labelT{countT}) = Dist3QcBTU.(labelP{countP}).(labelT{countT})*Price.Fuel*1.13; %$/hr
%             Dist3HXDOpCost.(labelP{countP}).(labelT{countT}) = Dist3HXDOpCost.(labelP{countP}).(labelT{countT})*(24*7*50)/1e6; %MM$/yr
%             
%             Dist3QrBTU.(labelP{countP}).(labelT{countT}) = Dist3Qr.(labelP{countP}).(labelT{countT})*0.947817/1e6; %MM Btu/hr
%             Dist3HXBOpCost.(labelP{countP}).(labelT{countT}) = Dist3QrBTU.(labelP{countP}).(labelT{countT})*Price.Fuel*0.93; %$/hr
%             Dist3HXBOpCost.(labelP{countP}).(labelT{countT}) = Dist3HXBOpCost.(labelP{countP}).(labelT{countT})*(24*7*50)/1e6; %MM$/yr
%             
%             
            %Sum of all three columns
            DistSteamEnergy.(labelP{countP}).(labelT{countT}) = ... %kJ/hr
                Dist1Qr.(labelP{countP}).(labelT{countT}) + ...
                Dist2Qr.(labelP{countP}).(labelT{countT});
            DistCost.(labelP{countP}).(labelT{countT}) = ...  $
                Dist1HXDCost.(labelP{countP}).(labelT{countT}) + ...
                Dist1HXBCost.(labelP{countP}).(labelT{countT}) + ...
                Dist2HXDCost.(labelP{countP}).(labelT{countT}) + ...
                Dist2HXBCost.(labelP{countP}).(labelT{countT}) + ...
                Dist1ColumnCost.(labelP{countP}).(labelT{countT}) + ...
                Dist2ColumnCost.(labelP{countP}).(labelT{countT});
            DistCost.(labelP{countP}).(labelT{countT}) = DistCost.(labelP{countP}).(labelT{countT})/1e6; %MM$
            DistOpCost.(labelP{countP}).(labelT{countT}) = ...
                Dist1HXDOpCost.(labelP{countP}).(labelT{countT}) + ...
                Dist2HXDOpCost.(labelP{countP}).(labelT{countT});
                Dist1HXBOpCost.(labelP{countP}).(labelT{countT}) + ...
                Dist2HXBOpCost.(labelP{countP}).(labelT{countT});
        end
        DistCost.(labelP{countP}).(labelT{countT}) = DistCost.(labelP{countP}).(labelT{countT})';
        DistOpCost.(labelP{countP}).(labelT{countT}) = DistOpCost.(labelP{countP}).(labelT{countT})';
        DistSteamEnergy.(labelP{countP}).(labelT{countT}) = DistSteamEnergy.(labelP{countP}).(labelT{countT})';
%         DistCost.(labelP{countP}).(labelT{countT}) = zeros(length(A),1); %ignoring cost
%         DistOpCost.(labelP{countP}).(labelT{countT}) = zeros(length(A),1); %ignoring operating cost
%         DistSteamEnergy.(labelP{countP}).(labelT{countT}) = zeros(length(A),1); %ignoring cost %kJ/hr
        
        %turbine
        gamma = 0.23;
        Pin = 1*2116.22; %lbf/ft^2
        Pout = P*2116.22; %lbf/ft^2
        Qin = (1/60)*0.0353147.*PTotal.(labelP{countP}).(labelT{countT}).*R*300/1; %ft^3/min
        hp = (3.03e-5/gamma)*Pin.*Qin*((Pout/Pin)^gamma - 1);
        bhp = hp/0.8;
        Fd = 1.15; %reciprocal, turbine
        Fc = Fd;
        TurbCost.(labelP{countP}).(labelT{countT}) = (MS/280)*517.5*(bhp.^0.82)*(2.11 + Fc);
        TurbCost.(labelP{countP}).(labelT{countT}) = TurbCost.(labelP{countP}).(labelT{countT})/1e6;
        TurbCost.(labelP{countP}).(labelT{countT}) = zeros(length(A),1);
        %turbine operating costs
        TurbOPCost.(labelP{countP}).(labelT{countT}) = -(bhp/0.9)*(1/1.341)*Price.Electricity*8400/1e6; % MM$/year
        TurbOPCost.(labelP{countP}).(labelT{countT}) = zeros(length(A),1); %ignoring operating cost due to gained electricity powering pumps
        
        %carbon operating costs
        ElectricityCarbon.(labelP{countP}).(labelT{countT}) = zeros(length(A),1);
        SteamEnergy.(labelP{countP}).(labelT{countT}) = HX1Q.(labelP{countP}).(labelT{countT}) + DistSteamEnergy.(labelP{countP}).(labelT{countT}); %kJ/hr
        FlowMethaneForSteam.(labelP{countP}).(labelT{countT}) = SteamEnergy.(labelP{countP}).(labelT{countT})/MethaneEnergy; %gmol/hr
        SteamCarbon.(labelP{countP}).(labelT{countT}) = MethaneCarbon*FlowMethaneForSteam.(labelP{countP}).(labelT{countT}); %gmol/hr
        SteamCarbon.(labelP{countP}).(labelT{countT}) = SteamCarbon.(labelP{countP}).(labelT{countT})*44.01/1000; %kg/hr
        SteamCarbon.(labelP{countP}).(labelT{countT}) = SteamCarbon.(labelP{countP}).(labelT{countT})*(24*50*7)/1000; %MT/year
        Carbon.(labelP{countP}).(labelT{countT}) = FurnaceCarbon.(labelP{countP}).(labelT{countT}) + ...
            ElectricityCarbon.(labelP{countP}).(labelT{countT}) + ...
            SteamCarbon.(labelP{countP}).(labelT{countT}); %including carbon costs, $
        CarbonCost.(labelP{countP}).(labelT{countT}) = Carbon.(labelP{countP}).(labelT{countT})*Price.CO2/1e6; %MM$
        %             CarbonCost.(label{count}) = zeros(length(A),1); %ignoring carbon costs
        
        % total upfront equipment costs
        Equip.(labelP{countP}).(labelT{countT}) = HXcost.(labelP{countP}).(labelT{countT}) + ...
            ReactorCost.(labelP{countP}).(labelT{countT}) + DistCost.(labelP{countP}).(labelT{countT}) + ...
            TurbCost.(labelP{countP}).(labelT{countT}) + FurnaceCost.(labelP{countP}).(labelT{countT}); %MM$
        FCI.(labelP{countP}).(labelT{countT}) = 2.28*Equip.(labelP{countP}).(labelT{countT});%MM$
        FCI.(labelP{countP}).(labelT{countT}) = FCImultiplier*FCI.(labelP{countP}).(labelT{countT}); %MM$ (for sensitivityanalysis)
        SU.(labelP{countP}).(labelT{countT}) = 0.1*FCI.(labelP{countP}).(labelT{countT});%MM$
        SU.(labelP{countP}).(labelT{countT}) = SUmultiplier*SU.(labelP{countP}).(labelT{countT}); %MM$ (for sensitivityanalysis)
        WC.(labelP{countP}).(labelT{countT}) = 2/12*FreshFeedM.(labelP{countP}).(labelT{countT})*Price.Ethanol/1000*24*7*50/1e6;%MM$/year
        WC.(labelP{countP}).(labelT{countT}) = WCmultiplier*WC.(labelP{countP}).(labelT{countT}); %MM$ (for sensitivityanalysis)
        TI.(labelP{countP}).(labelT{countT}) = FCI.(labelP{countP}).(labelT{countT}) + SU.(labelP{countP}).(labelT{countT}) + WC.(labelP{countP}).(labelT{countT});%MM$
        
        %loan amount
        loan.(labelP{countP}).(labelT{countT}) = TI.(labelP{countP}).(labelT{countT});
        %construction rate for first 2 years, 2% higher than corporate rate
        ConstructionInterest.(labelP{countP}).(labelT{countT}) = loan.(labelP{countP}).(labelT{countT})*(ConstructionInterestRate);
        
        %corporate bond cash flows, 3% (DOW/BASF), total amount of construction
        %costs and interest up to that point, 10 year bond life (simple
        %interest)
        bond.(labelP{countP}).(labelT{countT}) = loan.(labelP{countP}).(labelT{countT}) + ConstructionInterest.(labelP{countP}).(labelT{countT}) +...
            ConstructionInterest.(labelP{countP}).(labelT{countT})*(1.05);
        CorporateInterest.(labelP{countP}).(labelT{countT}) = bond.(labelP{countP}).(labelT{countT})*(0.03);
        
        %chemical revenue, operating costs, and total revenuebefore taxes
        Revenue.(labelP{countP}).(labelT{countT}) = (24*7*50)*(PMEtA.(labelP{countP}).(labelT{countT}).*Price.EtA + ...
            PMEthanol.(labelP{countP}).(labelT{countT}).*Price.Ethanol- FMEthanol.(labelP{countP}).(labelT{countT}).*Price.Ethanol)/(1e9); %MM$/yr
        RevenueProducts.(labelP{countP}).(labelT{countT}) = (24*7*50)*(PMEtA.(labelP{countP}).(labelT{countT}).*Price.EtA)/1e9; %MM$/yr
        AGS.(labelP{countP}).(labelT{countT}) = 0.1*RevenueProducts.(labelP{countP}).(labelT{countT});
        OperateCost.(labelP{countP}).(labelT{countT}) = DistOpCost.(labelP{countP}).(labelT{countT}) + ...
            HeatOPCost.(labelP{countP}).(labelT{countT}) + TurbOPCost.(labelP{countP}).(labelT{countT}) + ...
            CarbonCost.(labelP{countP}).(labelT{countT});
        % cost of manufacturing
        COM.(labelP{countP}).(labelT{countT}) = OperateCost.(labelP{countP}).(labelT{countT}) + AGS.(labelP{countP}).(labelT{countT});
        
        %salvage value, assuming 1% of initial equipment costs
        salvage.(labelP{countP}).(labelT{countT}) = Equip.(labelP{countP}).(labelT{countT})*0.01/2;
        
        %depreciation, linear over 10 years
        deprec.(labelP{countP}).(labelT{countT}) = (FCI.(labelP{countP}).(labelT{countT}) + SU.(labelP{countP}).(labelT{countT}) ...
            - salvage.(labelP{countP}).(labelT{countT}))./10; %/year
        tax = IncomeTaxRate;
        
        %%calculate cash flows
        %assign initial cash flows, assuming 60/40 split during construction
        cashflow.year1.(labelP{countP}).(labelT{countT}) = zeros(length(A),1) + loan.(labelP{countP}).(labelT{countT}); %land
        cashflow.year2.(labelP{countP}).(labelT{countT}) = -(0.6)*FCI.(labelP{countP}).(labelT{countT});
        %startup capital
        cashflow.year3.(labelP{countP}).(labelT{countT}) = -(0.4)*FCI.(labelP{countP}).(labelT{countT}) - SU.(labelP{countP}).(labelT{countT}) - ...
            WC.(labelP{countP}).(labelT{countT}) - ConstructionInterest.(labelP{countP}).(labelT{countT});
        %cashflows during  10 year depreciation and corporate bond period
        for year = 4:12
            yearlist = ['year' num2str(year)];
            profit.(yearlist).(labelP{countP}).(labelT{countT}) = (Revenue.(labelP{countP}).(labelT{countT}) - ...
                COM.(labelP{countP}).(labelT{countT}) - CorporateInterest.(labelP{countP}).(labelT{countT}) - ...
                deprec.(labelP{countP}).(labelT{countT}));
            cashflow.(yearlist).(labelP{countP}).(labelT{countT}) = (Revenue.(labelP{countP}).(labelT{countT}) - ...
                COM.(labelP{countP}).(labelT{countT}) - CorporateInterest.(labelP{countP}).(labelT{countT}) - ...
                deprec.(labelP{countP}).(labelT{countT})).*(1-tax) + ...
                deprec.(labelP{countP}).(labelT{countT});
        end
        cashflow.year13.(labelP{countP}).(labelT{countT}) = Revenue.(labelP{countP}).(labelT{countT}) - ...
            COM.(labelP{countP}).(labelT{countT}) - CorporateInterest.(labelP{countP}).(labelT{countT}) - ...
            loan.(labelP{countP}).(labelT{countT});
        %after depreciation and bond life
        for year = 14:18
            yearlist = ['year' num2str(year)];
            cashflow.(yearlist).(labelP{countP}).(labelT{countT}) = (Revenue.(labelP{countP}).(labelT{countT}) - ...
                COM.(labelP{countP}).(labelT{countT})).*(1-tax);
        end
        % Before tax cash flows
        for year = 4:13
            yearlist = ['year' num2str(year)];
            cashflowBT.(yearlist).(labelP{countP}).(labelT{countT}) = (Revenue.(labelP{countP}).(labelT{countT}) - ...
                COM.(labelP{countP}).(labelT{countT}) - CorporateInterest.(labelP{countP}).(labelT{countT}));
        end
        for year = 14:18
            yearlist = ['year' num2str(year)];
            cashflowBT.(yearlist).(labelP{countP}).(labelT{countT}) = (Revenue.(labelP{countP}).(labelT{countT}) - ...
                COM.(labelP{countP}).(labelT{countT}));
        end
        
        %salvage costs
        cashflow.year18.(labelP{countP}).(labelT{countT}) = cashflow.year18.(labelP{countP}).(labelT{countT}) + ...
            salvage.(labelP{countP}).(labelT{countT}) + WC.(labelP{countP}).(labelT{countT});
        %catalyst costs every 2 years
        CatalystVol.(labelP{countP}).(labelT{countT}) = Vreal.(labelP{countP}).(labelT{countT})*1000; %mL
        CatalystMass.(labelP{countP}).(labelT{countT}) = CatalystVol.(labelP{countP}).(labelT{countT})*0.5*6.31/1e6; %MT
        CatalystCost.(labelP{countP}).(labelT{countT}) = CatalystMass.(labelP{countP}).(labelT{countT})*Price.Cat/1e6; %MM$
        for year = 3:1:18
            yearlist = ['year' num2str(year)];
            cashflow.(yearlist).(labelP{countP}).(labelT{countT}) = cashflow.(yearlist).(labelP{countP}).(labelT{countT}) - ...
                CatalystCost.(labelP{countP}).(labelT{countT})*2;
        end
        cashflownoloans = cashflow;
        cashflownoloans.year1.(labelP{countP}).(labelT{countT}) = cashflownoloans.year1.(labelP{countP}).(labelT{countT}) - loan.(labelP{countP}).(labelT{countT});
        cashflownoloans.year13.(labelP{countP}).(labelT{countT}) = cashflownoloans.year13.(labelP{countP}).(labelT{countT}) + loan.(labelP{countP}).(labelT{countT});
        %%NPV calculation
        NPV.(labelP{countP}).(labelT{countT}) = zeros(length(A),1);
        for year = 1:18
            yearlist = ['year' num2str(year)];
            discount.(yearlist) = (1/((1+EnterpriseRate)^year));
            NPVcashflow.(yearlist).(labelP{countP}).(labelT{countT}) = cashflow.(yearlist).(labelP{countP}).(labelT{countT}).*discount.(yearlist);
            NPVcashflownoloan.(yearlist).(labelP{countP}).(labelT{countT}) = cashflownoloans.(yearlist).(labelP{countP}).(labelT{countT}).*discount.(yearlist);
            NPV.(labelP{countP}).(labelT{countT}) = NPV.(labelP{countP}).(labelT{countT}) + NPVcashflow.(yearlist).(labelP{countP}).(labelT{countT});
        end
        % cumulative cash flows
        CumulativeDiscountedCashFlowLoan.year0.(labelP{countP}).(labelT{countT}) = zeros(length(A),1);
        CumulativeDiscountedCashFlowNoLoan.year0.(labelP{countP}).(labelT{countT}) = zeros(length(A),1);
        CumulativeNonDiscountedCashFlowLoan.year0.(labelP{countP}).(labelT{countT}) = zeros(length(A),1);
        CumulativeNonDiscountedCashFlowNoLoan.year0.(labelP{countP}).(labelT{countT}) = zeros(length(A),1);
        
        for year = 1:18
            yearlist = ['year' num2str(year)];
            yearlistprevious = ['year' num2str(year-1)];
            CumulativeDiscountedCashFlowLoan.(yearlist).(labelP{countP}).(labelT{countT}) = ...
                CumulativeDiscountedCashFlowLoan.(yearlistprevious).(labelP{countP}).(labelT{countT})...
                + NPVcashflow.(yearlist).(labelP{countP}).(labelT{countT});
            CumulativeDiscountedCashFlowNoLoan.(yearlist).(labelP{countP}).(labelT{countT}) = ...
                CumulativeDiscountedCashFlowNoLoan.(yearlistprevious).(labelP{countP}).(labelT{countT})...
                + NPVcashflownoloan.(yearlist).(labelP{countP}).(labelT{countT});
            CumulativeNonDiscountedCashFlowLoan.(yearlist).(labelP{countP}).(labelT{countT}) = ...
                CumulativeNonDiscountedCashFlowLoan.(yearlistprevious).(labelP{countP}).(labelT{countT})...
                + cashflow.(yearlist).(labelP{countP}).(labelT{countT});
            CumulativeNonDiscountedCashFlowNoLoan.(yearlist).(labelP{countP}).(labelT{countT}) = ...
                CumulativeNonDiscountedCashFlowLoan.(yearlistprevious).(labelP{countP}).(labelT{countT})...
                + cashflownoloans.(yearlist).(labelP{countP}).(labelT{countT});
        end
        
        %%NPV %
        NPVpercent.(labelP{countP}).(labelT{countT}) = 100*NPV.(labelP{countP}).(labelT{countT})./(TI.(labelP{countP}).(labelT{countT})*15);
        
        %%ROI calculation
        totalprofit.(labelP{countP}).(labelT{countT}) = zeros(length(A),1);
        for year = 4:18
            yearlist = ['year' num2str(year)];
            totalprofit.(labelP{countP}).(labelT{countT}) = totalprofit.(labelP{countP}).(labelT{countT}) +...
                cashflowBT.(yearlist).(labelP{countP}).(labelT{countT});
        end
        avgprofit.(labelP{countP}).(labelT{countT}) = totalprofit.(labelP{countP}).(labelT{countT})./15;
        ROIBT.(labelP{countP}).(labelT{countT}) = 100.*avgprofit.(labelP{countP}).(labelT{countT})./TI.(labelP{countP}).(labelT{countT});
        %%IRR calculation
        %comment out fsolve for speed
        options = optimoptions(@fsolve,'DiffMaxChange',0.005,'FunctionTolerance',0.01,'Display','off');
        rate = zeros(length(NPV.(labelP{countP}).(labelT{countT})),1);
        for z = 2:length(NPV.(labelP{countP}).(labelT{countT}))
            if NPV.(labelP{countP}).(labelT{countT})(z) < 2
                rate0 = -3;
            else
                rate0 = 0.1;
            end
            p.countT = countT;
            p.countP = countP;
            p.labelT = labelT;
            p.labelP = labelP;
            p.z = z;
            p.cashflow = cashflownoloans;
            % use IRR calculation
                                rate(z) = fsolve(@(rate) IRRsolver(rate,p),rate0,options);
            % ignore IRR to speed up code
%             rate(z) = 0;
            IRR.(labelP{countP}).(labelT{countT})(z) = rate(z);
        end
    end
end

% for j = 10:5:30
%     P = j; %bar
%     countP = (j-10)/5 + 1;
%     for i = 300:20:380
%         countT = (i - 300)/20 + 1;
%         [maxes(countT),I1(countT)] = max(NPV.(labelP{countP}).(labelT{countT}));
%         maxes2(countT) = max(NPVpercent.(labelP{countP}).(labelT{countT}));
%         maxes3(countT) = max(IRR.(labelP{countP}).(labelT{countT}));
%     end
%     [maxNPV.(labelP{countP}),I] = max(maxes);
%     maxNPVpercent.(labelP{countP}) = max(maxes2);
%     maxIRR.(labelP{countP}) = max(maxes3);
%     % index.(labelP{countP}) = I1(I)
% end

% %cash flow diagram
% for year = 1:18
%     conversionindex = 25;
%     countT = 7;
%     yearlist = ['year' num2str(year)];
%     yearlistprevious = ['year' num2str(year-1)];
%     CumDiscCashFlow(year) = CumulativeDiscountedCashFlowLoan.(yearlist).(labelP{countP}).(labelT{countT})(conversionindex);
%     CumCashFlow(year) = CumulativeNonDiscountedCashFlowLoan.(yearlist).(labelP{countP}).(labelT{countT})(conversionindex);
%     CumDiscCashFlowNoLoan(year) = CumulativeDiscountedCashFlowNoLoan.(yearlist).(labelP{countP}).(labelT{countT})(conversionindex);
%     CumCashFlowNoLoan(year) = CumulativeNonDiscountedCashFlowNoLoan.(yearlist).(labelP{countP}).(labelT{countT})(conversionindex);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------plotting-------------------------%
countP = 1;
for z=1
    %     figure1 = figure('Color',[1 1 1]);
    %     figure2 = figure('Color',[1 1 1]);
    %     figure3 = figure('Color',[1 1 1]);
    figure4 = figure('Color',[1 1 1]);
    %     figure5 = figure('Color',[1 1 1]);
    %     figure6 = figure('Color',[1 1 1]);
    %     figure7 = figure('Color',[1 1 1]);
    %     figure8 = figure('Color',[1 1 1]);
    figure9 = figure('Color',[1 1 1]);
    %     figure10 = figure('Color',[1 1 1]);
    figure11 = figure('Color',[1 1 1]);
    figure12 = figure('Color',[1 1 1]);
    figure13 = figure('Color',[1 1 1]);
    figure14 = figure('Color',[1 1 1]);
    figure15 = figure('Color',[1 1 1]);
    figure16 = figure('Color',[1 1 1]);
    %     figure17 = figure('Color',[1 1 1]);
    
    %     figure(figure1)
    %     for i = 225:5:240
    %         countT = (i - 225)/5 + 1;
    %         hold on
    %         % plot mole fraction as function of tau
    %         plot(tau.(labelP{countP}).(labelT{countT}).*60,Ethanol.(labelP{countP}).(labelT{countT}),'Color',cc(countT,:),'LineWidth',1)
    %         xlabel('Residence Time [min]');ylabel('Molar Flow Rate Ethanol [gmol/hr]')
    %         %         axis([0 500 4 10])
    %     end
    %     axes1 = get(figure1,'CurrentAxes');
    %     hold(axes1,'on');
    %     box(axes1,'on');
    %     set(axes1,'FontWeight','bold');
    %     %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %     %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','NE','FontSize',11)
    %
    %     figure(figure2)
    %     for i = 225:5:240
    %         countT = (i - 225)/5 + 1;
    %         %plot molar flow rate PEtA vs residence time
    %         hold on
    %         plot((tau.(labelP{countP}).(labelT{countT})).*60,EtA.(labelP{countP}).(labelT{countT}),'Color',cc(countT,:),'LineWidth',1)
    %         xlabel('Residence Time [min]');ylabel('Molar Flow Rate Ethyl Acetate [gmol/hr]')
    %     end
    %     %     axis([0 500 0 2.2])
    %     axes1 = get(figure2,'CurrentAxes');
    %     hold(axes1,'on');
    %     box(axes1,'on');
    %     set(axes1,'FontWeight','bold');
    %     %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %     %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','best','FontSize',11)
    %
    %     figure(figure3)
    %     for i = 225:5:240
    %         countT = (i - 225)/5 + 1;
    %         %
    %         hold on
    %         plot(xt.(labelP{countP}).(labelT{countT}),sEtA.(labelP{countP}).(labelT{countT}),'Color',cc(countT,:),'LineWidth',1)
    %         xlabel('Ethanol Conversion');
    %         ylabel('Ethyl Acetate Selectivity')
    %     end
    %     %     axis([0 0.6 0.35 0.5])
    %     axes1 = get(figure3,'CurrentAxes');
    %     hold(axes1,'on');
    %     box(axes1,'on');
    %     set(axes1,'FontWeight','bold');
    %     %     yline = [-30:300]';
    %     %     oppoint = 0.25;
    %     %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %     %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','NE','FontSize',11)
    %
    %     %
    figure(figure4)
    for i = 225:5:240
        countT = (i - 225)/5 + 1;
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),Revenue.(labelP{countP}).(labelT{countT}),'Color',cc(countT,:),'Marker','sq')
        xlabel('Ethanol Conversion');ylabel('Revenue [MM$]')
    end
    %     axis([0 0.6 25 40])
    axes1 = get(figure4,'CurrentAxes');
    hold(axes1,'on');
    box(axes1,'on');
    set(axes1,'FontWeight','bold');
    %     yline = [-30:300]';
    %     oppoint = 0.25;
    %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','Best','FontSize',11)
    
    %
    %     figure(figure5)
    %     for i = 225:5:240
    %         countT = (i - 225)/5 + 1;
    %         hold on
    %         plot(xt.(labelP{countP}).(labelT{countT}),Vreal.(labelP{countP}).(labelT{countT})./1000,'Color',cc(countT,:),'LineWidth',1)
    %         xlabel('Ethanol Conversion');
    %         ylabel('Reactor Volume [m^3]')
    %     end
    %     %     axis([0 0.6 0 3*10^3])
    %     axes1 = get(figure5,'CurrentAxes');
    %     hold(axes1,'on');
    %     box(axes1,'on');
    %     set(axes1,'FontWeight','bold');
    %     %     yline = [-30:30000]';
    %     %     oppoint = 0.25;
    %     %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %     %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','NE','FontSize',11)
    %
    %     %
    %     figure(figure7)
    %     for i = 225:5:240
    %         countT = (i - 225)/5 + 1;
    %         hold on
    %         plot(xt.(labelP{countP}).(labelT{countT}),(FEthanol.(labelP{countP}).(labelT{countT})-PEthanol.(labelP{countP}).(labelT{countT})),'Color',cc(countT,:),'LineWidth',1)
    %         xlabel('Ethanol Conversion');ylabel('Fresh Feed Ethanol [gmol/hr]')
    %     end
    %     %     axis([0 0.6 160000 240000]);
    %     axes1 = get(figure7,'CurrentAxes');
    %     hold(axes1,'on');
    %     box(axes1,'on');
    %     set(axes1,'FontWeight','bold');
    %     %     yline = [-30:300]';
    %     %     oppoint = 0.25;
    %     %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %     %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','Best','FontSize',11)
    %
    %     %
    figure(figure9)
    for i = 240
        countT = (i - 225)/5 + 1;
        yyaxis left
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),NPV.(labelP{countP}).(labelT{countT}),'Color','k','LineStyle','none','Marker','sq','LineWidth',1.5)
        xlabel('Ethanol Conversion')
        ylabel('NPV_{project} [MM$]')
                axis([0.2 0.6 61 65.5])
        
        yyaxis right
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),NPVpercent.(labelP{countP}).(labelT{countT}),'Color','k','LineStyle','none','Marker','o','LineWidth',1.5)
        ylabel('NPV_%')
                axis([0.2 0.6 5 20])
    end
    axes1 = get(figure9,'CurrentAxes');
    hold(axes1,'on');
    box(axes1,'on');
    set(axes1,'FontWeight','bold','YColor','k');
    yyaxis left
    axes1 = get(figure9,'CurrentAxes');
    set(axes1,'FontWeight','bold','YColor','k');
        yline = [-300:300]';
        oppoint = 0.45;
        plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     legend(['340' char(176) 'C'],...
    %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','NE','FontSize',11)

    %     %
    %     figure(figure10)
    %     for i = 240
    %         countT = (i - 225)/5 + 1;
    %         hold on
    %         plot(xt.(labelP{countP}).(labelT{countT}),PEther.(labelP{countP}).(labelT{countT})./PTotal.(labelP{countP}).(labelT{countT}),...
    %             xt.(labelP{countP}).(labelT{countT}),PEtA*ones(length(xt.(labelP{countP}).(labelT{countT})),1)./PTotal.(labelP{countP}).(labelT{countT}),...
    %             xt.(labelP{countP}).(labelT{countT}),PEthanol.(labelP{countP}).(labelT{countT})./PTotal.(labelP{countP}).(labelT{countT}),'LineWidth',1)
    %         yline = [-30:300]';
    %         oppoint = 0.25;
    %         plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %         legend('x_{Ether}','x_{Ethyl Acetate}','x_{Ethanol}','Location','best','FontSize',11)
    %         xlabel('Ethanol Conversion');
    %         ylabel('Reactor Effluent Mole Fraction')
    %     end
    %     %     axis([0 0.6 0 1])
    %     axes1 = get(figure10,'CurrentAxes');
    %     hold(axes1,'on');
    %     box(axes1,'on');
    %     set(axes1,'FontWeight','bold');
    
    %
    figure(figure11)
    for i = 225:5:240
        countT = (i - 225)/5 + 1;
        hold on
        plot(NPV.(labelP{countP}).(labelT{countT}),NPVpercent.(labelP{countP}).(labelT{countT}),'Color',cc(countT,:),'marker','sq')
    end
    %     axis([40 80 0 30])
    ylabel('NPV_%')
    xlabel('NPV_{project} [MM$]')
    axes1 = get(figure11,'CurrentAxes');
    hold(axes1,'on');
    box(axes1,'on');
    set(axes1,'FontWeight','bold');
    %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','Best','FontSize',11)
    
    %
    figure(figure12)
    for i = 240
        countT = (i - 225)/5 + 1;
        hold on
        IRRplot.(labelP{countP}).(labelT{countT}) = 100*IRR.(labelP{countP}).(labelT{countT});
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),IRRplot.(labelP{countP}).(labelT{countT}),'Color','k','Marker','sq','LineWidth',1.5,'LineStyle','none')
    end
    %     axis([0 0.6 0 45])
    ylabel('IRR [%]')
    xlabel('Ethanol Conversion')
    axes1 = get(figure12,'CurrentAxes');
    hold(axes1,'on');
    box(axes1,'on');
    set(axes1,'FontWeight','bold');
    %     yline = [-30:300]';
    %     oppoint = 0.25;
    %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','NE','FontSize',11)
    
    %
    figure(figure13)
    for i = 225:5:240
        countT = (i - 225)/5 + 1;
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),ROIBT.(labelP{countP}).(labelT{countT}),'Color',cc(countT,:),'Marker','sq')
    end
    %     axis([0 0.6 -10 70])
    ylabel('ROI_{BT} [%]')
    xlabel('Ethanol Conversion')
    axes1 = get(figure13,'CurrentAxes');
    hold(axes1,'on');
    box(axes1,'on');
    set(axes1,'FontWeight','bold');
    %     yline = [-30:300]';
    %     oppoint = 0.25;
    %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     legend(['300' char(176) 'C'],['320' char(176) 'C'],['340' char(176) 'C'],...
    %         ['360' char(176) 'C'],['380' char(176) 'C'],'Location','NE','FontSize',11)
    
    %
    figure(figure14)
    for i = 225:5:240
        countT = (i - 225)/5 + 1;
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),TI.(labelP{countP}).(labelT{countT}),'Color',cc(countT,:),'Marker','sq')
    end
    %     axis([0 0.6 0 100])
    ylabel('TCI [MM$]')
    xlabel('Ethanol Conversion')
    axes1 = get(figure14,'CurrentAxes');
    hold(axes1,'on');
    box(axes1,'on');
    set(axes1,'FontWeight','bold');
    %     yline = [-30:300]';
    %     oppoint = 0.25;
    %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    
    %
    figure(figure15)
    for i = 240
        countT = (i - 225)/5 + 1;
        yyaxis left
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),ROIBT.(labelP{countP}).(labelT{countT}),'Color','k','LineStyle',':','Marker','sq','LineWidth',1.5)
        xlabel('Ethanol Conversion')
        ylabel('ROI_{BT} [%]')
                axis([0.2 0.6 20 50])
        
        yyaxis right
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),TI.(labelP{countP}).(labelT{countT}),'Color','k','LineStyle',':','Marker','o','LineWidth',1.5)
        ylabel('TCI [MM $]')
                axis([0.2 0.6 35 70])
    end
    axes1 = get(figure15,'CurrentAxes');
    hold(axes1,'on');
    box(axes1,'on');
    set(axes1,'FontWeight','bold','YColor','k');
    yyaxis left
    axes1 = get(figure15,'CurrentAxes');
    set(axes1,'FontWeight','bold','YColor','k');
        legend(['ROI_{BT}'],['TI'],['Operating Point'],'Location','NW','FontSize',11)
yyaxis right
        yline = [-30:300]';
        oppoint = 0.45;
        plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
            legend(['ROI_{BT}'],['TI'],['Operating Point'],'Location','NW','FontSize',11)

    %     %
    %     figure(figure6)
    %     for i = 240
    %         countT = (i - 225)/5 + 1;
    %         hold on
    %         plot(xt.(labelP{countP}).(labelT{countT}),PEther.(labelP{countP}).(labelT{countT})./1000,...
    %             xt.(labelP{countP}).(labelT{countT}),PEtA*ones(1,length(xt.(labelP{countP}).(labelT{countT})))./1000,...
    %             xt.(labelP{countP}).(labelT{countT}),PEthanol.(labelP{countP}).(labelT{countT})./1000,'LineWidth',1.5)
    %         hold on
    %         plot(xt.(labelP{countP}).(labelT{countT}),FEthanol.(labelP{countP}).(labelT{countT})./1000,'LineWidth',1.5)
    %         plot(xt.(labelP{countP}).(labelT{countT}),FreshFeed.(labelP{countP}).(labelT{countT})./1000,'LineWidth',1.5)
    %         %     yline = [-30:3000]';
    %         %     oppoint = 0.25;
    %         %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %         legend('P_{Ether}','P_{Ethyl Acetate}','R_{Ethanol}','F_{Ethanol}','F_{Fresh Ethanol}','Location','best','FontSize',11)
    %         xlabel('Toluene Conversion');ylabel('Molar Flowrate [kgmol/hr]')
    %         %     axis([0 0.6 0 10*max(PEtA)./1000])
    %         axes1 = get(figure6,'CurrentAxes');
    %         hold(axes1,'on');
    %         box(axes1,'on');
    %         set(axes1,'FontWeight','bold');
    %     end
    %
    %     figure(figure8)
    %     %only one temp as PTotal is independent of temp for a given conversion
    %     for i = 240
    %         countT = (i - 225)/5 + 1;
    %         hold on
    %         plot(xt.(labelP{countP}).(labelT{countT}),PTotal.(labelP{countP}).(labelT{countT}),'Color',cc(countT,:),'LineWidth',1)
    %         xlabel('Ethanol Conversion');ylabel('Total Flow to Sep Unit [gmol/hr]')
    %         %     axis([0 0.6 0 1e7])
    %         axes1 = get(figure8,'CurrentAxes');
    %         hold(axes1,'on');
    %         box(axes1,'on');
    %         set(axes1,'FontWeight','bold');
    %         %     yline = [-30:30000000]';
    %         %     oppoint = 0.25;
    %         %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     end
    %
    
    %
    figure(figure16)
    countT = 4;
    for j = 2:2:10
        P = j; %bar
        countP = (j-2)/2 + 1;
        yyaxis left
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),NPV.(labelP{countP}).(labelT{countT}),'Color',cc(countP,:),'LineStyle','-','Marker','sq','LineWidth',1)
        xlabel('Toluene Conversion')
        ylabel('NPV_{project} [MM$]')
        %         axis([0 0.6 -10 80])
        
        yyaxis right
        hold on
        plot(xt.(labelP{countP}).(labelT{countT}),NPVpercent.(labelP{countP}).(labelT{countT}),'Color',cc(countP,:),'LineStyle',':','Marker','o','LineWidth',1)
        ylabel('NPV_%')
        %         axis([0 0.6 0 30])
    end
    axes1 = get(figure16,'CurrentAxes');
    hold(axes1,'on');
    box(axes1,'on');
    set(axes1,'FontWeight','bold','YColor','k');
    yyaxis left
    axes1 = get(figure16,'CurrentAxes');
    set(axes1,'FontWeight','bold','YColor','k');
    %     yline = [-300:300]';
    %     oppoint = 0.25;
    %         plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %
    %     legend(['10 bar'],['20 bar'],...
    %         ['30 bar'],'Location','NE','FontSize',11)
    %
    %
    %         %
    %     figure(figure17)
    %     countT = 5;
    %     for i = 1:2:5
    %         countP = i;
    %         hold on
    %         plot(xt.(labelP{countP}).(labelT{countT}),Vreal.(labelP{countP}).(labelT{countT})./1000,'Color',cc(countP,:),'LineWidth',1)
    %         xlabel('Toluene Conversion');
    %         ylabel('Reactor Volume [m^3]')
    %     end
    %     axis([0 0.6 0 1*10^2])
    %     axes1 = get(figure17,'CurrentAxes');
    %     hold(axes1,'on');
    %     box(axes1,'on');
    %     set(axes1,'FontWeight','bold');
    %     yline = [-30:30000]';
    %     oppoint = 0.25;
    %     plot(zeros(length(yline),1)+oppoint,yline,'LineStyle','--','Color','k','Linewidth',1,'Marker','none')
    %     legend(['10 bar'],['20 bar'],['30 bar'],'Location','NE','FontSize',11)
    
    
    
end
% -------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% toc

function F = IRRsolver(rate,p)
countP = p.countP;
countT = p.countT;
labelP = p.labelP;
labelT = p.labelT;
z = p.z;
cashflow = p.cashflow;
NPVIRR.(labelP{countP}).(labelT{countT})(z) = 0;
for year = 1:18
    yearlist = ['year' num2str(year)];
    discount.(yearlist) = (1/((1 + rate)^year));
    NPVIRRcashflow.(yearlist).(labelP{countP}).(labelT{countT}) = cashflow.(yearlist).(labelP{countP}).(labelT{countT})(z).*discount.(yearlist);
    NPVIRR.(labelP{countP}).(labelT{countT})(z) = NPVIRR.(labelP{countP}).(labelT{countT})(z) + NPVIRRcashflow.(yearlist).(labelP{countP}).(labelT{countT});
end
F(1) = NPVIRR.(labelP{countP}).(labelT{countT})(z);
end