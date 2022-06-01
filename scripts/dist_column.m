function [Qc,Qr,Height,Diameter,Nreal,r,s,VoverF] = dist_column(F,zF,rmin,P,T,lambdaD,lambdaB,MvFeed,MvDist,MvBottoms,aij)
R = 8.3145e-5;
%mole fractions are for the light component
%sat liquid feed
q = 1;

D = zF.*F; %assumes zF is mole fraction of LK AND all lighter components
B = (1-zF).*F;

r = 1.5.*rmin;
s = (D./B).*(r + q) - (1 - q); %reboiler ratio
% i = LK, j = HK
fiD = 0.997;
fjB = 0.997;
fjD = 1 - fiD;
fiB = 1 - fjB;

Nmin = log((fiD./fiB).*(fjB./fjD))/log(aij); %Fenske equation
xi = 0.75.*(1-((r-rmin)./(r+1)).^0.5688); %dummy variable
Ntheoretical = (xi + Nmin)./(1-xi); %FUG method, p 136
Nreal = 2.*Ntheoretical; %(can replace with OConnell's correlation if desired, p 260)

Ht = 60; %cm
Hmin = 3.*Ht; %cm
Height = Hmin + Ht.*Nreal; %cm
Height = Height./100; %m

Vt = (r+1).*D;
Vmin = (rmin+1).*D;
Vb = s.*B;
V = Vb;
VoverF = V/F;

rhol = 870; %kg/m^3

for i = 1:length(MvFeed)
Mv(i,1) = (MvFeed(i)+MvDist(i)+MvBottoms(i))/3;
end

% ideal gas law
rhogfeed = P.*MvFeed./(R.*T);
rhogdist = P.*MvDist./(R.*T);
rhogbottoms = P.*MvBottoms./(R.*T);
for i = 1:length(rhogfeed)
rhog(i,1) = (rhogfeed(i)+rhogdist(i)+rhogbottoms(i))/3;
end
c0 = 439;
phi = 0.6;

Area = (Mv./sqrt(rhol.*rhog))*(1./(phi.*c0)).*(.1./0.8).*V;
Diameter = 2.*sqrt(Area./pi);

Qc = lambdaD.*Vt; %kJ/hr
Qr = lambdaB.*Vb; %kJ/hr
end
