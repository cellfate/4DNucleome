%% 取得间隔更密  为了画分界线 Figure S17 & Figure 3
rates = 0.005:0.001:0.3;
rhoBS_temp = zeros(size(rates,2),size(rates,2),size(rates,2));
rhoBF_temp = zeros(size(rates,2),size(rates,2),size(rates,2));

for idx = 1:size(rates,2)
    for jdx = 1:size(rates,2)
        params.simulated_on      = true;
        params.EP_flag           = true;
        params.compute_TR        = false;
        params.attraction_coef   = 0.1; % K, E-P communication coefficient
        params.enhancer_index    = 25; % enhancer index number
        params.promoter_index    = 75; % promoter index number
        
        % Downstream: Transcriptional bursting
        params.k_on1_max         = rates(idx); % burst initiation rate
        params.k_on1             = 0.002; % burst initiation rate
        params.k_recruitment_max = rates(jdx); % burst polymerase recruitment rate
        params.k_recruitment     = 0.009; % burst polymerase recruitment rate
        params.k_release_max     = rates; % burst polymerase pause release rate
        params.k_release         = 0.008; % burst polymerase pause release rate
        params.k_off2            = 0.006;% off1 to off2
        params.k_on2             = 0.002;% off2 to off1
        params.k_off1            = 0.009;% rec to off1
        params.k_off3            = 0.002;% rec to off2
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
        rhoBS_temp(idx,jdx,:) = BS_max./BS_min;
        rhoBF_temp(idx,jdx,:) = CT_min./CT_max;
    end
end


params.diffusion_const = 4e-4;
params.friction_coef = 10;
params.spring_const = 1;
params.promoter_index = 75;
params.enhancer_index = 25;

params.distance_T     = 0.1; 
params.distance_05    = 0.2;
params.H              = 3;
%% 确定速率
struct_rates = [];
slopeBSBF_V2 = zeros(15*5);
for i = 1:1:15*5
    for j = 1:1:15*5
        A = ((rhoBF_temp > i/5-0.01 & rhoBF_temp < i/5+0.01) & (rhoBS_temp > j/5-0.01 & rhoBS_temp < j/5+0.01));
        tot_num = sum(sum(sum(A)));
        if tot_num == 0
            A = ((rhoBF_temp > i/5-0.05 & rhoBF_temp < i/5+0.05) & (rhoBS_temp > j/5-0.05 & rhoBS_temp < j/5+0.05));
            tot_num = sum(sum(sum(A)));
        end
        if tot_num == 0
            A = ((rhoBF_temp > i/5-0.1 & rhoBF_temp < i/5+0.1) & (rhoBS_temp > j/5-0.1 & rhoBS_temp < j/5+0.1));
            tot_num = sum(sum(sum(A)));
        end        
        if tot_num == 0
            A = ((rhoBF_temp > i/5-0.15 & rhoBF_temp < i/5+0.15) & (rhoBS_temp > j/5-0.15 & rhoBS_temp < j/5+0.15));
            tot_num = sum(sum(sum(A)));
        end
        if tot_num ~= 0
            [x,~,~] = find(A==1);
            count = 1;
            for idx = 1:size(x)
                struct_rates.rhoBSBF(i,j).num(count).on1 = rates(x(idx));
                [y,z] = find(reshape(A(x(idx),:,:),[size(rates,2),size(rates,2)]) == 1);
                struct_rates.rhoBSBF(i,j).num(count).rec = rates(y(1));
                struct_rates.rhoBSBF(i,j).num(count).rel = rates(z(1));
                A(x(idx),y(1),z(1)) = 0;
                struct_rates.rhoBSBF(i,j).num(count).bs = rhoBS_temp(x(idx),y(1),z(1));
                struct_rates.rhoBSBF(i,j).num(count).bf = rhoBF_temp(x(idx),y(1),z(1));
                
                
                params.k_on1_max = rates(x(idx));
                params.k_recruitment_max = rates(y(1));
                params.k_release_max = rates(z(1));
                %% 计算斜率
                params.attraction_coef = 1;
                Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
                    (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
                dt = 0.01;
                d_EP = 0.01:dt:5;
                P = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
                
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
                % BSV1
                burst_size_theorV1 = (recruitment.*release./(recruitment.*off3+(release+off3).*off1));
                % ONPDFV1
                total_on_time_theorV1 = (recruitment+release+off3)/(recruitment.*off3+(release+off3).*off1);
                % OFFV1
                total_off_time_theorV1 = (off2+on2)./(on2.*on1);
                % BSV2
                burst_size_theorV2 = sum((lambda_recruitment.*lambda_release./(lambda_recruitment.*lambda_off3+...
                    (lambda_release+lambda_off3).*lambda_off1)).*P.*dt);
                % ONV2
                total_on_time_theorV2 = sum((lambda_recruitment+lambda_release+lambda_off3)./...
                    (lambda_recruitment.*lambda_off3+(lambda_release+lambda_off3).*lambda_off1).*P.*dt);
                % OFFV2
                total_off_time_theorV2 = sum((lambda_off2+lambda_on2)./(lambda_on2.*lambda_on1).*P.*dt);
                
                % mu
                minlambda = min([params.k_on1,params.k_recruitment,params.k_release]);
                maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
                    params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
                    99,1, 'first'))/params.friction_coef;
                mu = minlambda/maxvelocity;
                burst_size_theor_ep1 = 1./(1+mu).*burst_size_theorV1+mu./(1+mu).*burst_size_theorV2;
                total_on_time_theor = 1./(1+mu).*total_on_time_theorV1+mu./(1+mu).*total_on_time_theorV2;
                total_off_time_theor = 1./(1+mu).*total_off_time_theorV1+mu./(1+mu).*total_off_time_theorV2;
                burst_frequency_theor_ep1 = 1./(total_on_time_theor+total_off_time_theor);
                
                %% 1. kep = 0.001;
                params.attraction_coef = 0.001;
                Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
                    (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
                dt = 0.01;
                d_EP = 0.01:dt:5;
                P = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
                
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
                
                % V2
                burst_size_theorV2 = sum((lambda_recruitment.*lambda_release./...
                    (lambda_recruitment.*lambda_off3+...
                    (lambda_release+lambda_off3).*lambda_off1)).*P.*dt);
                % ON
                total_on_time_theorV2 = sum((lambda_recruitment+lambda_release+lambda_off3)./...
                    (lambda_recruitment.*lambda_off3+(lambda_release+lambda_off3).*lambda_off1).*P.*dt);
                % OFF
                total_off_time_theorV2 = sum((lambda_off2+lambda_on2)./(lambda_on2.*lambda_on1).*P.*dt);
                
                % mu
                minlambda = min([params.k_on1,params.k_recruitment,params.k_release]);
                maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
                    params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
                    99,1, 'first'))/params.friction_coef;
                mu = minlambda/maxvelocity;
                burst_size_theor_ep002 = 1./(1+mu).*burst_size_theorV1+mu./(1+mu).*burst_size_theorV2;
                total_on_time_theor = 1./(1+mu).*total_on_time_theorV1+mu./(1+mu).*total_on_time_theorV2;
                total_off_time_theor = 1./(1+mu).*total_off_time_theorV1+mu./(1+mu).*total_off_time_theorV2;
                burst_frequency_theor_ep002 = 1./(total_on_time_theor+total_off_time_theor);
                
                dlogBS_dlogk = (log10(burst_size_theor_ep1) - log10(burst_size_theor_ep002));
                dlogBF_dlogk = (log10(burst_frequency_theor_ep1) - log10(burst_frequency_theor_ep002));
                struct_rates.rhoBSBF(i,j).num(count).slopeBSBF_V1 = dlogBS_dlogk/dlogBF_dlogk;
                count = count + 1;
                
            end
            struct_rates.rhoBSBF(i,j).tot_count = count - 1;
        end
        slopeBSBF_V2(i,j) = log10(j/(1+(j-1)*0.0224))/log10(i+(1-i)*0.0224);
    end
end


%%
slopeBSBF_V1_mean = zeros(15*5);
slopeBSBF_V1_max = zeros(15*5);
slopeBSBF_V1_min = zeros(15*5);
max_rates = [];
for i = 1:1:15*5
    for j = 1:1:15*5
        if struct_rates.rhoBSBF(i,j).tot_count ~= 0
            temp = [];
            for idx = 1:struct_rates.rhoBSBF(i,j).tot_count
                temp(end+1) = struct_rates.rhoBSBF(i,j).num(idx).slopeBSBF_V1;
            end
            [~,maxindex] = max(temp);
            [~,minindex] = min(temp);
            slopeBSBF_V1_max(i,j) = max(temp);
            slopeBSBF_V1_min(i,j) = min(temp);
            slopeBSBF_V1_mean(i,j) = mean(temp);
            if 2*slopeBSBF_V1_min(i,j) < slopeBSBF_V1_max(i,j)
                max_rates(end+1,1) = i;
                max_rates(end,2) = j;
                max_rates(end,3) = struct_rates.rhoBSBF(i,j).num(maxindex).on1;
                max_rates(end,4) = struct_rates.rhoBSBF(i,j).num(maxindex).rec;
                max_rates(end,5) = struct_rates.rhoBSBF(i,j).num(maxindex).rel;
                max_rates(end,6) = struct_rates.rhoBSBF(i,j).num(maxindex).slopeBSBF_V1;
                max_rates(end,7) = struct_rates.rhoBSBF(i,j).num(minindex).on1;
                max_rates(end,8) = struct_rates.rhoBSBF(i,j).num(minindex).rec;
                max_rates(end,9) = struct_rates.rhoBSBF(i,j).num(minindex).rel;
                max_rates(end,10) = struct_rates.rhoBSBF(i,j).num(minindex).slopeBSBF_V1;
                max_rates(end,11) = slopeBSBF_V2(i,j);
            end
        end
    end
end

%%

slopeBSBF_V2 = zeros(1400);
m = 1; n = 1;
for i = 1.01:0.01:15
    for j = 1.01:0.01:15
        slopeBSBF_V2(m,n) = log10(j/(1+(j-1)*0.0224))/log10(i+(1-i)*0.0224);
        n = n + 1;
    end
    m = m + 1; n = 1;
end
figure
imagesc(slopeBSBF_V2(101:end,1:end)>1)
axis xy
axis square
set(gca,'TickLength',[0.02,0.025]);

X = 10:75;
Y = zeros(1,size(X,2));
j = 1;
for i = 10:75
    Y(j) = find(slopeBSBF_V1_max(i,:)>1,1,"first") - 1;
    if j >=24 && j<=52 %max
        Y(j) = Y(j)+2;
    end
    j = j+ 1;
end

% mean
x = 10:75;
p1 =   7.055e-07 ;
       p2 =  -0.0001602  ;
       p3 =      0.0134 ;
       p4 =     -0.5051  ;
       p5 =       8.816  ;
       p6 =      -40.27 ;
Z = p1.*x.^5 + p2.*x.^4 + p3.*x.^3 + p4.*x.^2 + p5.*x + p6;
figure
plot(x,Z)
axis([10 75 5 75])
axis square
hold on
% min
x = 10:75;
p1 =    7.14e-07  ;
       p2 =   -0.000158 ;
       p3 =     0.01256  ;
       p4 =     -0.4438  ;
       p5 =       7.679 ;
       p6 =      -34.09 ;
Z = p1.*x.^5 + p2.*x.^4 + p3.*x.^3 + p4.*x.^2+p5.*x.^1 + p6;
plot(x,Z)
axis([10 75 5 75])
axis square

% max
x = 10:75;
p1 =   6.973e-07 ;
       p2 =  -0.0001523 ;
       p3 =     0.01249 ;
       p4 =     -0.4697  ;
       p5 =        8.17 ;
       p6 =      -36.58 ;
Z = p1.*x.^5 + p2.*x.^4 + p3.*x.^3 + p4.*x.^2+p5.*x.^1 + p6;
plot(x,Z)
axis([10 75 5 75])
axis square
set(gca,'Layer','top','TickLength',[0.02 0.025],'XTick',[25 50 75],...
    'XTickLabel',{'5','10','15'},'YTick',[25 50 75],'YTickLabel',{'5','10','15'});