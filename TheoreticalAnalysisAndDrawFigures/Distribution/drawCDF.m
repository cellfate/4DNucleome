addpath(genpath('.\ResultsKEP'));
%% distributions of KEP Figure S13A-D
K_coefficient =  [ -1.8000   -1.6000   -1.4000   -1.2000   -1.0000  -0.8000   -0.6000   -0.4000   -0.2000   -0.0000 ];
figure1 = figure;
for i = 1:size(K_coefficient,2)
    filename = sprintf('ResultsKEP\\K_%d.mat',10^(K_coefficient(i)));
    load(filename);
    hold on
    plot(results.TheroPDF.t,cumsum(results.TheroPDF.PCTPDF)./100);
end
box on
set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);

figure1 = figure;
for i = 1:size(K_coefficient,2)
    filename = sprintf('ResultsKEP\\K_%d.mat',10^(K_coefficient(i)));
    load(filename);
    hold on
    plot(results.TheroPDF.BS_Bin,cumsum(results.TheroPDF.BSPDF));
end
box on
set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);

figure1 = figure;
for i = 1:size(K_coefficient,2)
    filename = sprintf('ResultsKEP\\K_%d.mat',10^(K_coefficient(i)));
    load(filename);
    hold on
    plot(results.TheroPDF.t,cumsum(results.TheroPDF.PONPDF)./100);
end
box on
set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);

%% dG
addpath(genpath('.\ResultsdG'));
enhancer_idx = [51,55,60,65,70,75,80,85,90,95];
promoter_idx = [49,45,40,35,30,25,20,15,10,5];
figure1 = figure;
for i = 1:size(promoter_idx,2)
    filename = sprintf('ResultsdG\\E_%d_P_%d',enhancer_idx(i),promoter_idx(i));
    load(filename);
    hold on
    plot(results.TheroPDF.t,cumsum(results.TheroPDF.PCTPDF)./100);
end
box on
set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);

figure1 = figure;
for i = 1:size(promoter_idx,2)
    filename = sprintf('ResultsdG\\E_%d_P_%d',enhancer_idx(i),promoter_idx(i));
    load(filename);
    hold on
    plot(results.TheroPDF.BS_Bin,cumsum(results.TheroPDF.BSPDF));
end
box on
set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);

figure1 = figure;
for i = 1:size(promoter_idx,2)
    filename = sprintf('ResultsdG\\E_%d_P_%d',enhancer_idx(i),promoter_idx(i));
    load(filename);
    hold on
    plot(results.TheroPDF.t,cumsum(results.TheroPDF.PONPDF)./100);
end
box on
set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);