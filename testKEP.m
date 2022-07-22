% clear;clc;
% close all;
%% current parallel pool
if isempty(gcp('nocreate'))
    parpool(24);
end
addpath(genpath('.\SimulationAndStatisticalAnalysis'));
addpath(genpath('.\TheoreticalAnalysisAndDrawFigures'));
K_coefficient = -1.8:0.2:0;
for i = 1:1:size(K_coefficient,2)
    %% Parameters setting
    % Upstream: Chromatin conformation
    input_options.simulated_on      = true;
    input_options.EP_flag           = true;
    input_options.compute_TR        = false;
    input_options.attraction_coef   = 10^(K_coefficient(i)); % K, E-P communication coefficient
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
    
    input_options.simulation_num    = 100; % simulation number
    input_options.simulation_reaction_step = 1000000; %
    input_options.simulation_time   = 2000000; % [s]
    
    input_options.result_base_folder = fullfile(pwd, 'ResultsKEP');
    input_options.filename = sprintf('//K_%d',input_options.attraction_coef);
    if exist([input_options.result_base_folder,input_options.filename],'file') == 0
        mkdir(input_options.result_base_folder,input_options.filename);
    end
    
    %% Parallel computing
    tic;
    if input_options.simulated_on
        parfor s_idx = 1:input_options.simulation_num
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
    results.TheroPDF = AnalyseBurstPDF(params);
    filename = sprintf('//K_%d.mat',input_options.attraction_coef);
    save([input_options.result_base_folder,filename],'results');
    timerVal = toc;
    X = ['Analysing time:',num2str(timerVal)];
    disp(X)
end







