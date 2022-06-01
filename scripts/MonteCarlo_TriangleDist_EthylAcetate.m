clear
close all
N = 10000;
tic
NPVlist = zeros(1,N);
for k = 1:N
%     tic
    %range of effects [-%, +%]
    EnterpriseRate = 0.1;
    GeneralRange = [-.1 .1];
    
    EthRange = [-0.01 0.01];
    DERange = [-.25 .5];
    EARange = [-.32 .24];
    FuelRange = GeneralRange;
    CO2Range = GeneralRange;
    CatalystRange = GeneralRange;
    InterestRateRange = [-.5 .5];
    IncomeTaxRateRange = [-.5 .5];
    FCIRange = [-.1 .25];
    SURange = [-.1 1.0];
    WCRange = [-.2 .5];
    
    Price.Eth = 500;%[$/MT]
    Price.DE = 800;%[$/MT]
    Price.EA = 1200;%[$/MT]
    Price.Fuel = 2.5;%[$/MM BTU]
    Price.CO2 = 40; %$/MT
    Price.Cat = 10000; %$/MT
    InterestRate = 0.05;
    IncomeTaxRate = 0.27;
    FCImultiplier = 1;
    SUmultiplier = 1;
    WCmultiplier = 1;
    
    M(1,1) = Price.Eth;
    M(1,2) = Price.DE;
    M(1,3) = Price.EA;
    M(1,4) = Price.Fuel;
    M(1,5) = Price.CO2;
    M(1,6) = Price.Cat;
    M(1,7) = InterestRate;
    M(1,8) = IncomeTaxRate;
    M(1,9) = FCImultiplier;
    M(1,10) = SUmultiplier;
    M(1,11) = WCmultiplier;
    
    M([2 3],1) = EthRange;
    M([2 3],2) = DERange;
    M([2 3],3) = EARange;
    M([2 3],4) = FuelRange;
    M([2 3],5) = CO2Range;
    M([2 3],6) = CatalystRange;
    M([2 3],7) = InterestRateRange;
    M([2 3],8) = IncomeTaxRateRange;
    M([2 3],9) = FCIRange;
    M([2 3],10) = SURange;
    M([2 3],11) = WCRange;
    L = length(M(1,:));
    %assume triangle distributions for all ranges
    for q = 1:L
        % lower bound
        a(q) = (1+M(2,q))*M(1,q);
        % base value
        b(q) = M(1,q);
        % upper bound
        c(q) = (1+M(3,q))*M(1,q);
        pd = makedist('Triangular','a',a(q),'b',b(q),'c',c(q));
        M(4,q) = random(pd,1);
    end
    
    Price.Eth = M(4,1);%[$/MT]
    Price.DE = M(4,2);%[$/MT]
    Price.EA = M(4,3);%[$/MT]
    Price.Fuel = M(4,4);%[$/MM BTU]
    Price.CO2 = M(4,5); %$/MT
    Price.Cat = M(4,6); %$/MT
    InterestRate = M(4,7);
    IncomeTaxRate = M(4,8);
    FCImultiplier = M(4,9);
    SUmultiplier = M(4,10);
    WCmultiplier = M(4,11);
        
    NPV = Sensitivity_Hysys_EthylAcetate(Price,EnterpriseRate,InterestRate,...
        IncomeTaxRate,FCImultiplier,SUmultiplier,WCmultiplier);
    NPVlist(k) = NPV;
end

% make cdf of NPVlist
orderedNPVlist = sort(NPVlist);
prob = zeros(length(orderedNPVlist),1);
for z = 1:length(orderedNPVlist)
    prob(z) = (z-1)/length(orderedNPVlist);
end


%plot histogram
figure1 = figure('Color',[1 1 1]);
histogram(NPVlist,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',.3)
ylabel('\boldmath{$Probability$}','Interpreter','latex','FontSize',12);
xlabel('\boldmath{$NPV, [MM\$]$}','Interpreter','latex','FontSize',12);
axes1 = get(figure1,'CurrentAxes');
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontWeight','bold');

%plot cdf
figure1 = figure('Color',[1 1 1]);
plot(orderedNPVlist,prob,'Marker','s','Color','k','LineStyle','none','LineWidth',1)
ylabel('\boldmath{$Probability$}','Interpreter','latex','FontSize',12);
xlabel('\boldmath{$NPV_{project}$}','Interpreter','latex','FontSize',12);
axes1 = get(figure1,'CurrentAxes');
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontWeight','bold');
% toc