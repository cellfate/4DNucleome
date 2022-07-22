% clear;clc;
% close all;
%% current parallel pool
if isempty(gcp('nocreate'))
    parpool(24);
end
addpath(genpath('.\SimulationAndStatisticalAnalysis'));
addpath(genpath('.\TheoreticalAnalysisAndDrawFigures'));
enhancer_idx = [51,55,60,65,70,80,85,90,95];
promoter_idx = [49,45,40,35,30,20,15,10,5];
for i = 1:size(enhancer_idx,2)
simulated_on = true;
%% Parameters setting
% Upstream: Chromatin conformation
input_options.simulated_on      = true;
input_options.EP_flag           = true;
input_options.compute_TR        = true;
input_options.attraction_coef   = 0.1; % K, E-P communication coefficient
input_options.enhancer_index    = enhancer_idx(i); % enhancer index number
input_options.promoter_index    = promoter_idx(i); % promoter index number

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

input_options.simulation_num    = 1000; % simulation number
input_options.simulation_reaction_step = 500000; %
input_options.simulation_time   = 4000000; % [s]

input_options.result_base_folder = fullfile(pwd, 'ResultsdG');
input_options.filename = sprintf('//E_%d_P_%d',input_options.enhancer_index,input_options.promoter_index);
if exist([input_options.result_base_folder,input_options.filename],'file') == 0
    mkdir(input_options.result_base_folder,input_options.filename);
end

%% Parallel computing
tic;
if input_options.simulated_on
    for s_idx = 1:input_options.simulation_num
        % Load parameters
        params = ParametersBurst(input_options);
        % simulation
        result = SimulateBurst(params,s_idx);
    end
end
timerVal = toc;
X = ['Total simulation time:',num2str(timerVal)];
disp(X)
%% Analysing Burst
% Only parameters are needed, and data are loaded from the .mat files
tic;
params = ParametersBurst(input_options);
results = AnalyseBurst(params);
% results.TheroPDF = AnalyseBurstPDF(params);
filename = sprintf('//E_%d_P_%d',input_options.enhancer_index,input_options.promoter_index);
save([input_options.result_base_folder,filename],'results');
timerVal = toc;
X = ['Analysing time:',num2str(timerVal)];
disp(X)
end







