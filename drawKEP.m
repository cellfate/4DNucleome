% change kep Figure 2
K_coefficient =  [ -1.8000   -1.6000   -1.4000   -1.2000   -1.0000  -0.8000   -0.6000   -0.4000   -0.2000   -0.0000 ];
output_size = size(K_coefficient,2);

burst_size = zeros(1,output_size);
burst_frequency = zeros(1,output_size);
burst_size_theor = zeros(1,output_size);
burst_frequency_theor = zeros(1,output_size);
addpath(genpath('.\ResultsKEP'));
%% loading data
for i = 1:1:output_size

    filename = sprintf('ResultsKEP\\K_%d.mat',10^(K_coefficient(i)));
    load(filename);
    burst_size(1,i) = results.burst_size;
    burst_frequency(1,i) = results.burst_frequency;
    burst_size_theor(1,i) =  results.burst_size_theor;
    burst_frequency_theor(1,i) = results.burst_frequency_theor;
    
end
%%
params = results.simulate(1).data.result.params;
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


dlogBS_dlogk = (log10(1.9991) - log10(0.797844521461438))./(log10(1)-log10(0.02));
dlogBF_dlogk = (log10(0.0026) - log10(7.924586845114539e-04))./(log10(1)-log10(0.02));

bs = log10(BS_max) + dlogBS_dlogk.*(K_coefficient);
bf = -log10(TOFF_max+TON_max) + dlogBF_dlogk.*(K_coefficient);
%% figure  % kep变化使用
%% burst size
figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
plot(K_coefficient,log10(burst_size_theor));
hold on
scatter((K_coefficient),log10(burst_size));
 plot(K_coefficient, bs);
box on
% axis([-1.7 0.1 -0.2 0.4])
set(gca,'TickLength',[0.02,0.025]);



%% burst frequency
figure2 = figure;
set(figure2,'position',[300 400 280 190],'Name','bf');
hold on
plot((K_coefficient),log10(burst_frequency_theor));
scatter((K_coefficient),log10(burst_frequency));
 plot(K_coefficient, bf);
% axis([-1.7 0.1 -3.2 -2.5])
box on
set(gca,'TickLength',[0.02,0.025]);
