% distributions and properties of E-P spatial distance
%% Figure S10C left
F = figure;
for i = 0.1:0.1:1
    Kaa = sqrt(0.004*...
        (1/(abs(50))...
        +i)^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:0.8;
    EPPDF_Thero = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    mbcdf = cumsum(EPPDF_Thero)/100;
    hold on
    set(F,'position',[300 400 280 190]);
    plot((d_EP),(mbcdf));
    set(gca,'TickLength',[0.02,0.025]);
    box on
end
%% Figure S10C left
F = figure;
for i = 1:10:90
    Kaa = sqrt(0.004*...
        (1/(abs(i))...
        +0.1)^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:0.8;
    EPPDF_Thero = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    mbcdf = cumsum(EPPDF_Thero)/100;
    hold on
    set(F,'position',[300 400 280 190]);
    plot((d_EP),(mbcdf));
    set(gca,'TickLength',[0.02,0.025]);
    box on
end

%% Figure S10E left mean vs contact prob. power law
F = figure;
kep = 0.01:0.01:1;
ds_mean = zeros(1,length(kep));
connet_prob = zeros(1,length(kep));
for i = 1:length(kep)
    Kaa = sqrt(0.004*...
        (1/(abs(50))...
        +kep(i))^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:5;
    EP_Bin_Thero = d_EP;
    EPPDF_Thero = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    mbcdf = cumsum(EPPDF_Thero)/100;
    ds_mean(i) = 2*Kaa*sqrt(2/pi);
    connet_prob(i) = mbcdf(10);
end
A = ds_mean;
B = connet_prob;

dG = 2:1:99;
ds_mean = zeros(1,length(dG));
connet_prob = zeros(1,length(dG));
for i = 1:length(dG)
    Kaa = sqrt(0.004*...
        (1/(abs(dG(i)))...
        +0.1)^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:5;
    EP_Bin_Thero = d_EP;
    EPPDF_Thero = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    mbcdf = cumsum(EPPDF_Thero)/100;
    ds_mean(i) = 2*Kaa*sqrt(2/pi);
    connet_prob(i) = mbcdf(10);
end
A = [A,ds_mean];
B = [B,connet_prob];
set(F,'position',[300 400 280 190]);
scatter(log10(A(1:5:end)),log10(B(1:5:end)),30,'filled')
x = 0.1:0.001:0.6;
y = 0.003404.*x.^(-2.253);
hold on
plot(log10(x),log10(y))
set(gca,'TickLength',[0.02,0.025]);
box on


%% Figure S10E right dg vs encouter
F = figure;
dG = 1:1:99;
connet_prob = zeros(1,length(dG));
for i = 1:length(dG)
    Kaa = sqrt(0.004*...
        (1/(abs(dG(i)))...
        +0.1)^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:5;
    EP_Bin_Thero = d_EP;
    EPPDF_Thero = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    mbcdf = cumsum(EPPDF_Thero)/100;
    connet_prob(i) = mbcdf(10);
    
end
hold on
set(F,'position',[300 400 280 190]);
scatter(log10(dG(1:1:end)),log10(connet_prob(1:1:end)),30,'filled');
x = 0:1:90;
y = 0.5102.*x.^(-0.6796);

hold on
plot(log10(x),log10(y))
axis([0 2 -2 0])
set(gca,'TickLength',[0.02,0.025]);
box on


%% Figure S10F left
F = figure;
for i = 0.1:0.1:1
    Kaa = sqrt(0.004*...
        (1/(abs(50))...
        +i)^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:5;
    EPPDF_Thero = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    mbcdf = cumsum(EPPDF_Thero)/100;
    e_d = [2:1:200];
    for j = 1:size(e_d,2)
        encounter(j) = mbcdf(e_d(j))/((4/3)*pi*(d_EP(e_d(j))/2)^3);
    end
    hold on
    set(F,'position',[300 400 280 190]);
    plot(log10(d_EP(e_d)/2),log10(encounter));
    set(gca,'TickLength',[0.02,0.025]);
    box on
end

%% Figure S10F right
F = figure;
for i = 2:10:90
    Kaa = sqrt(0.004*...
        (1/(abs(i))...
        +0.1)^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:5;
    EPPDF_Thero = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    mbcdf = cumsum(EPPDF_Thero)/100;
    e_d = [2:1:200];
    for j = 1:size(e_d,2)
        encounter(j) = mbcdf(e_d(j))/((4/3)*pi*(d_EP(e_d(j))/2)^3);
    end
    hold on
    set(F,'position',[300 400 280 190]);
    plot(log10(d_EP(e_d)/2),log10(encounter));
    set(gca,'TickLength',[0.02,0.025]);
    box on
end

