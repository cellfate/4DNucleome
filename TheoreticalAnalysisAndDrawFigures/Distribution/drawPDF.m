addpath(genpath('.\ResultsKEP'));
%% pdf for a certain parameter for Figure S13 E-H
K_coefficient1 =  -1.6;
filename = sprintf('ResultsKEP\\K_%d.mat',10^(K_coefficient1));
load(filename);

figBS = figure;
set(figBS,'position',[300 400 200 190]);
% h = histogram(results.mRNA,max(results.mRNA),'Normalization','pdf');
hold on
bar(results.TheroPDF.BS_Bin,results.TheroPDF.BSPDF)
set(gca,'TickLength',[0.02,0.025]);
axis([-1 11 0 0.5])


figOFF = figure;
h = histogram(results.total_off_time,100,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.POFFPDF)

figON = figure;
h = histogram(results.total_on_time,100,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.PONPDF)

figCT = figure;
h = histogram(results.total_off_time+results.total_on_time,100,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.PCTPDF)




