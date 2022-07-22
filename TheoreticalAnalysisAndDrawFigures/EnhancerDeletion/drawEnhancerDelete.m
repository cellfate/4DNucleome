addpath(genpath('.\ResultsKEP'));
%% Figure S11B-E
figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
h = histogram(results.mRNA,max(results.mRNA),'Normalization','pdf');
x = h.BinEdges(1:end-1);
y = h.Values;
bar(x,y);
hold on
plot(results.TheroPDF.BS_Bin,results.TheroPDF.BSPDF)
box on
axis([-0.75 6.75 0 0.6])
set(gca,'TickLength',[0.02,0.025]);

figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
h = histogram(results.total_on_time,700,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.PONPDF)
box on
set(gca,'TickLength',[0.02,0.025]);

figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
h = histogram(results.total_off_time,700,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.POFFPDF)
box on
set(gca,'TickLength',[0.02,0.025]);

figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
h = histogram(results.total_ct_time,700,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.PCTPDF)
box on
set(gca,'TickLength',[0.02,0.025]);

%% Figure S11F
burst_size = [0.861388636592996 0.669838444687842,2.005136931836254,1.509254915340982,1.186126657660390,1.062882243622233];
burst_size_theor = [0.781766870680841 0.666666666666667,1.999141717873241,1.482492769558345,1.092421907641133,0.962847273723692];
burst_frequency = [7.079842337392362e-04 4.568121161596675e-04,0.002621485836021,0.001955470638435,0.001295588805264,0.001057131295274];
burst_frequency_theor = [6.852129714055832e-04 4.595744680851064e-04,0.002603674210711,0.001930349051850,0.001266718734590,0.001029619665199];

figure1 = figure;
x = [burst_size burst_size_theor];
y = [burst_frequency burst_frequency_theor];
x = log10(x);
y = log10(y);
hold on
scatter(x(1:12),y(1:12));
box on
set(gca,'TickLength',[0.02,0.025]);
t = -0.3:0.1:0.7;
s = 1.573845445616245.*t-3.0689;
hold on
plot(t,s)

t = -0.3:0.1:0.7;
s = 1.398667904896260.*t-3.0139;
hold on
plot(t,s)
axis square
set(figure1,'position',[300 400 280 190]);
axis([-0.3 0.7 -3.5 -2.5])












