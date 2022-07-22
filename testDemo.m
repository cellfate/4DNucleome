clear;clc;
close all;
%% current parallel pool
if isempty(gcp('nocreate'))
    parpool(4);
end
addpath(genpath('.\SimulationAndStatisticalAnalysis'));
addpath(genpath('.\TheoreticalAnalysisAndDrawFigures'));
%% Parameters setting
% Upstream: Chromatin conformation
input_options.simulated_on      = true;
input_options.EP_flag           = true;
input_options.compute_TR        = false;
input_options.attraction_coef   = 1; % K, E-P communication coefficient
input_options.enhancer_index    = 25; % enhancer index number
input_options.promoter_index    = 75; % promoter index number

% Downstream: Transcriptional bursting
input_options.k_on1_max         = 0.020; % burst initiation rate
input_options.k_on1             = 0.002; % burst initiation rate
input_options.k_recruitment_max = 0.024; % burst polymerase recruitment rate
input_options.k_recruitment     = 0.009; % burst polymerase recruitment rate
input_options.k_release_max     = 0.025; % burst polymerase pause release rate
input_options.k_release         = 0.008; % burst polymerase pause release rate
input_options.k_off2            = 0.006;% off1 to off2
input_options.k_on2             = 0.002;% off2 to off1
input_options.k_off1            = 0.009;% rec to off1
input_options.k_off3            = 0.002;% rec to off2

input_options.simulation_num    = 4; % simulation number
input_options.simulation_reaction_step = 1000000; %
input_options.simulation_time   = 200000; % [s]

input_options.result_base_folder = fullfile(pwd, 'Results');
input_options.filename = sprintf('//K_%d',input_options.attraction_coef);
if exist([input_options.result_base_folder,input_options.filename],'file') == 0
    mkdir(input_options.result_base_folder,input_options.filename);
end


%% Parallel computing
if input_options.simulated_on
    tic;
    parfor s_idx = 1:input_options.simulation_num
        % Load parameters
        params = ParametersBurst(input_options);
        % simulation
        result = SimulateBurst(params,s_idx);
    end
    timerVal = toc;
    X = ['Total simulation time:',num2str(timerVal)];
    disp(X)
end
%% Analysing Burst
% Only parameters are needed, and data are loaded from the .mat files
tic;
params = ParametersBurst(input_options);
results = AnalyseBurst(params);
results.TheroPDF = AnalyseBurstPDF(params);
filename = sprintf('//K_%d.mat',input_options.attraction_coef);
save([input_options.result_base_folder,filename],'results');
timerVal = toc;
X = ['Analysing time:',num2str(timerVal)];
disp(X)

%% Figure 
figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
h = histogram(results.mRNA,max(results.mRNA),'Normalization','pdf');
x = h.BinEdges(1:end-1);
y = h.Values;
bar(x,y);
hold on
plot(results.TheroPDF.BS_Bin,results.TheroPDF.BSPDF,'LineWidth',1)
box on
axis([-0.75 15 0 0.4])
set(gca,'TickLength',[0.02,0.025]);
ylabel('Prob.');
xlabel('m');

figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
h = histogram(results.total_on_time,200,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.PONPDF,'LineWidth',1)
box on
axis([-20 1000 0 0.01])
set(gca,'TickLength',[0.02,0.025]);
ylabel('PDF');
xlabel('t');

figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
h = histogram(results.total_off_time,400,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.POFFPDF,'LineWidth',1)
box on
axis([-20 1000 0 0.02])
set(gca,'TickLength',[0.02,0.025]);
ylabel('PDF');
xlabel('t');

figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
h = histogram(results.total_on_time + results.total_off_time,300,'Normalization','pdf');
hold on
plot(results.TheroPDF.t,results.TheroPDF.PCTPDF,'LineWidth',1)
box on
set(gca,'TickLength',[0.02,0.025]);
axis([-40 2000 0 0.003])
ylabel('PDF');
xlabel('t');





