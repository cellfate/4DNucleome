% change dG Figure 4
enhancer_idx = [51,55,60,65,70,75,80,85,90,95];
promoter_idx = [49,45,40,35,30,25,20,15,10,5];
output_size = size(enhancer_idx,2);

burst_size = zeros(1,output_size);
burst_frequency = zeros(1,output_size);
burst_size_theor = zeros(1,output_size);
burst_frequency_theor = zeros(1,output_size);
addpath(genpath('.\ResultsdG'));
%% loading data
for i = 1:1:output_size

    filename = sprintf('ResultsdG\\E_%d_P_%d', enhancer_idx(i),promoter_idx(i));
    load(filename);
    burst_size(1,i) = results.burst_size;
    burst_frequency(1,i) = results.burst_frequency;
    burst_size_theor(1,i) =  results.burst_size_theor;
    burst_frequency_theor(1,i) = results.burst_frequency_theor;
    
end

%% figure  非对数画图
%% burst size
figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
plot((enhancer_idx - promoter_idx),(burst_size_theor));
hold on
scatter((enhancer_idx - promoter_idx),(burst_size));
box on
axis([-5 100 0.8 2.4])
set(gca,'TickLength',[0.02,0.025]);

errorbar((enhancer_idx - promoter_idx),burst_size,burst_size_std*0.1,'LineWidth',1);


%% burst frequency
figure2 = figure;
set(figure2,'position',[300 400 280 190],'Name','bf');
hold on
plot((enhancer_idx - promoter_idx),(burst_frequency_theor));
scatter((enhancer_idx - promoter_idx),(burst_frequency));
box on
set(gca,'TickLength',[0.02,0.025]);
errorbar((enhancer_idx - promoter_idx),burst_frequency,burst_frequency_std*1,'LineWidth',1);


%% 对数线性性
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


dlogBS_dlogk = (log10(2.0111) - log10(1.0924))./(log10(1)-log10(50));
dlogBF_dlogk = (log10(0.0026) - log10(0.0012))./(log10(1)-log10(60));

bs = log10(BS_max) + dlogBS_dlogk.*(log10(abs(enhancer_idx - promoter_idx)));
bf = -log10(TOFF_max+TON_max) + dlogBF_dlogk.*(log10(abs(enhancer_idx - promoter_idx)));

%% burst size
figure1 = figure;
set(figure1,'position',[300 400 280 190],'Name','bs');
plot(log10(enhancer_idx - promoter_idx),log10(burst_size_theor));
hold on
scatter(log10(enhancer_idx - promoter_idx),log10(burst_size));
 plot(log10(enhancer_idx - promoter_idx), bs);
box on
 axis([0.2 2 -0 0.3])
set(gca,'TickLength',[0.02,0.025]);



%% burst frequency
figure2 = figure;
set(figure2,'position',[300 400 280 190],'Name','bf');
hold on
plot(log10(enhancer_idx - promoter_idx),log10(burst_frequency_theor));
scatter(log10(enhancer_idx - promoter_idx),log10(burst_frequency));
 plot(log10(enhancer_idx - promoter_idx), bf);
% axis([-1.7 0.1 -3.2 -2.5])
box on
set(gca,'TickLength',[0.02,0.025]);

