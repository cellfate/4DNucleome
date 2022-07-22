%% approximation and power-law for KEP Figure S14DE
K_coefficient =  -2:0.01:0;
output_size = size(K_coefficient,2);

burst_size = zeros(1,output_size);
burst_frequency = zeros(1,output_size);
burst_size_theor = zeros(1,output_size);
burst_frequency_theor = zeros(1,output_size);
burst_size_theor_a = zeros(1,output_size);
burst_frequency_theor_a = zeros(1,output_size);

addpath(genpath('.\Results'));
filename = sprintf('Results\\K_%d.mat',1.000000e-01);
load(filename);
params = results.simulate(1).data.result.params;

for i = 1: output_size
    params.attraction_coef = 10^(K_coefficient(i));
    
    Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
        (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:5;
    d_T = params.distance_T;
    d_05 = params.distance_05;
    P = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    F_EP_05 = erf(d_05/(sqrt(2)*Kaa)) - sqrt(2/pi)*(Kaa)^(-1)...
        *d_05*exp(-d_05^(2)/(2*(Kaa)^(2)));
    minlambda = min([params.k_on1,params.k_recruitment,params.k_release]);
    maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
        params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
        99,1, 'first'))/params.friction_coef;
    mu = minlambda/maxvelocity;
    %% theor
    
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
    
    burst_size_theor(1,i) = 1./(1+mu).*burst_size_theorV1+mu./(1+mu).*burst_size_theorV2;
    total_on_time_theor = 1./(1+mu).*total_on_time_theorV1+mu./(1+mu).*total_on_time_theorV2;
    total_off_time_theor = 1./(1+mu).*total_off_time_theorV1+mu./(1+mu).*total_off_time_theorV2;
    burst_frequency_theor(1,i) = 1./(total_on_time_theor+total_off_time_theor);
    
    %% approximation
    %% fast
    lambda_on1 = params.k_on1_max*F_EP_05 + params.k_on1*(1-F_EP_05);
    lambda_rec = params.k_recruitment_max*F_EP_05 + params.k_recruitment*(1-F_EP_05);
    lambda_rel = params.k_release_max*F_EP_05 + params.k_release*(1-F_EP_05);
    lambda_on2 = params.k_on2; lambda_off2 = params.k_off2;
    lambda_off1 = params.k_off1;lambda_off3 = params.k_off3;
    % bs
    burst_size_theorV1 = (lambda_rec.*lambda_rel./...
        (lambda_rec.*lambda_off3+...
        (lambda_rel+lambda_off3).*lambda_off1));
    % ON
    total_on_time_theorV1 = (lambda_rec+lambda_rel+lambda_off3)./...
        (lambda_rec.*lambda_off3+(lambda_rel+lambda_off3).*lambda_off1);
    % OFF
    total_off_time_theorV1 = (lambda_off2+lambda_on2)./(lambda_on2.*lambda_on1);
    total_ct_time_theorV1 = total_on_time_theorV1 + total_off_time_theorV1;
    %% slow
    BS_max = (params.k_recruitment_max.*params.k_release_max./...
        (params.k_recruitment_max.*params.k_off3+...
        (params.k_release_max+params.k_off3).*params.k_off1));
    BS_min = (params.k_recruitment.*params.k_release./...
        (params.k_recruitment.*params.k_off3+...
        (params.k_release+params.k_off3).*params.k_off1));
    TON_max = (params.k_recruitment_max+params.k_release_max+params.k_off3)./...
        (params.k_recruitment_max.*params.k_off3+(params.k_release_max+...
        params.k_off3).*params.k_off1);
    TON_min = (params.k_recruitment+params.k_release+params.k_off3)./...
        (params.k_recruitment.*params.k_off3+(params.k_release+...
        params.k_off3).*params.k_off1);
    TOFF_max = (1+params.k_off2./params.k_on2)./params.k_on1_max;
    TOFF_min = (1+params.k_off2./params.k_on2)./params.k_on1;
    CT_max = TON_max + TOFF_max;
    CT_min = TON_min + TOFF_min;
    % bs
    burst_size_theorV2 = BS_max*F_EP_05 + BS_min*(1-F_EP_05);
    % ct
    total_ct_time_theorV2 =  CT_max*F_EP_05 + CT_min*(1-F_EP_05);
    
    burst_size_theor_a(1,i) = 1./(1+mu).*burst_size_theorV1+mu./(1+mu).*burst_size_theorV2;
    total_ct_time_theor = 1./(1+mu).*total_ct_time_theorV1+mu./(1+mu).*total_ct_time_theorV2;
    burst_frequency_theor_a(1,i) = 1./total_ct_time_theor;
    
    %% power law
    dlogBS_dlogk = (log10(2.0051) - log10(0.7978))./(log10(1)-log10(0.02));
    dlogBF_dlogk = (log10(0.0026) - log10(7.1642e-04))./(log10(1)-log10(0.02));
    
    bs = log10(BS_max) + dlogBS_dlogk.*(K_coefficient);
    bf = -log10(TOFF_max+TON_max) + dlogBF_dlogk.*(K_coefficient);
end
%% figure  % kep变化使用
%% burst size
figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
plot(K_coefficient,log10(burst_size_theor));
hold on
scatter((K_coefficient),log10(burst_size));
plot(K_coefficient, log10(burst_size_theor_a));
plot(K_coefficient, bs);
box on
 axis([-2.05 0.05 -0.2 0.4])
set(gca,'TickLength',[0.02,0.025]);


%% burst frequency
figure2 = figure;
set(figure2,'position',[300 400 280 190],'Name','bf');
hold on
plot((K_coefficient),log10(burst_frequency_theor));
scatter((K_coefficient),log10(burst_frequency));
plot(K_coefficient,log10(burst_frequency_theor_a));
plot(K_coefficient, bf);
 axis([-2.05 0.05 -3.4 -2.5])
box on
set(gca,'TickLength',[0.02,0.025]);




