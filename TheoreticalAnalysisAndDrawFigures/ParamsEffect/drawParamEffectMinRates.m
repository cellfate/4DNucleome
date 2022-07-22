%% Effects of model parameters (min) on burst Figure S15D-F
K_coefficient =  -2:0.01:0;
on1_min =  linspace(0.0001,0.01,10);
output_size = size(K_coefficient,2);
output_size1 = size(on1_min,2);

burst_size_theor = zeros(output_size1,output_size);
burst_frequency_theor = zeros(output_size1,output_size);
addpath(genpath('.\ResultsKEP'));
filename = sprintf('ResultsKEP\\K_%d.mat',10^(0));
load(filename);
params = results.simulate(1).data.result.params;

for i = 1: output_size1
    input_options.EP_flag           = true;
    input_options.attraction_coef   = 0.1; % K, E-P communication coefficient
    input_options.enhancer_index    = 25; % enhancer index number
    input_options.promoter_index    = 75; % promoter index number
    
    % Downstream: Transcriptional bursting
    input_options.k_on1_max         = 0.020; % burst initiation rate
    input_options.k_on1             = 0.010; % burst initiation rate
    input_options.k_recruitment_max = 0.024; % burst polymerase recruitment rate
    input_options.k_recruitment     = 0.010; % burst polymerase recruitment rate
    input_options.k_release_max     = 0.025; % burst polymerase pause release rate
    input_options.k_release         = (on1_min(i)); % burst polymerase pause release rate
    input_options.k_off2            = 0.006;% off1 to off2
    input_options.k_on2             = 0.002;% off2 to off1
    input_options.k_off1            = 0.009;% rec to off1
    input_options.k_off3            = 0.002;% rec to off2
      params = ParametersBurst(input_options);
    %     params.k_on1_max = 10^on1_max(i);
    %     params.k_release_max = 10^rec_max(i);
%     params.k_off1 = 10^koff1(i);
    for j = 1:output_size
        params.attraction_coef = 10^(K_coefficient(j));
        Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
            (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
        dt = 0.01;
        d_EP = 0.01:dt:5;
        P = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
%         F_EP_05 = erf(d_05/(sqrt(2)*Kaa)) - sqrt(2/pi)*(Kaa)^(-1)...
%             *d_05*exp(-d_05^(2)/(2*(Kaa)^(2)));
        minlambda = min([params.k_on1,params.k_recruitment,params.k_release]);
        maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
            params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
            99,1, 'first'))/params.friction_coef;
        mu = minlambda/maxvelocity;
        %% theor
        d_T = params.distance_T;
        d_05 = params.distance_05;
        Hill = (d_EP > d_T).* (1./(1+((d_EP-d_T)./(d_05-d_T)).^params.H));
        Hill(1:d_T/dt) = 0;
        lambda_on1 = (d_EP <= d_T).*params.k_on1_max + (d_EP > d_T).*...
            (params.k_on1 + (params.k_on1_max -params.k_on1).*Hill);
        lambda_recruitment = (d_EP <= d_T).*params.k_recruitment_max + (d_EP > d_T).*...
            (params.k_recruitment + (params.k_recruitment_max -params.k_recruitment).*Hill);
        lambda_release = (d_EP <= d_T).*params.k_release_max + (d_EP > d_T).*...
            (params.k_release+ (params.k_release_max -params.k_release).*Hill);
        lambda_on2 = params.k_on2; lambda_off2 = params.k_off2;
        lambda_off1 = params.k_off1;lambda_off3 = params.k_off3;
        on1 = sum(lambda_on1.*P.*dt);
        recruitment = sum(lambda_recruitment.*P.*dt);
        release = sum(lambda_release.*P.*dt);
        on2 = sum(lambda_on2.*P.*dt);
        off2 = sum(lambda_off2.*P.*dt);
        off1 = sum(lambda_off1.*P.*dt);
        off3 = sum(lambda_off3.*P.*dt);
        % BS
        burst_size_theorV1 = (recruitment.*release./(recruitment.*off3+(release+off3).*off1));
        % ONPDF
        total_on_time_theorV1 = (recruitment+release+off3)/(recruitment.*off3+(release+off3).*off1);
        % OFF
        total_off_time_theorV1 = (off2+on2)./(on2.*on1);
        % tr
        traveling_ratio_theor = params.elongation_time*release;
        
        % bs
        burst_size_theorV2 = sum((lambda_recruitment.*lambda_release./...
            (lambda_recruitment.*lambda_off3+...
            (lambda_release+lambda_off3).*lambda_off1)).*P.*dt);
        % ON
        total_on_time_theorV2 = sum((lambda_recruitment+lambda_release+lambda_off3)./...
            (lambda_recruitment.*lambda_off3+(lambda_release+lambda_off3).*lambda_off1).*P.*dt);
        % OFF
        total_off_time_theorV2 = sum((lambda_off2+lambda_on2)./(lambda_on2.*lambda_on1).*P.*dt);
        
        burst_size_theor(i,j) = 1./(1+mu).*burst_size_theorV1+mu./(1+mu).*burst_size_theorV2;
        total_on_time_theor = 1./(1+mu).*total_on_time_theorV1+mu./(1+mu).*total_on_time_theorV2;
        total_off_time_theor = 1./(1+mu).*total_off_time_theorV1+mu./(1+mu).*total_off_time_theorV2;
        burst_frequency_theor(i,j) = 1./(total_on_time_theor+total_off_time_theor);
    end
    
end
%% figure  % kep变化使用
%% burst size
figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
plot(K_coefficient,log10(burst_size_theor));
hold on
box on
% axis([-1.7 0.1 -0.2 0.4])
set(gca,'TickLength',[0.02,0.025]);

% errorbar(K_coefficient,burst_size,burst_size_std*0.1,'LineWidth',1);


%% burst frequency
figure2 = figure;
set(figure2,'position',[300 400 280 190],'Name','bf');
hold on
plot((K_coefficient),log10(burst_frequency_theor));
% axis([-1.7 0.1 -3.2 -2.5])
box on
set(gca,'TickLength',[0.02,0.025]);
% errorbar(K_coefficient,burst_frequency_off,burst_frequency_off_std*0.3,'LineWidth',1);





