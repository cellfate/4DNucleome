%% logarithmic gains (two partial derivatives) for kep Figure S14F-H
clear;
dG = 10.^linspace(0,3,500); 
%  rates = 0.003:0.001:0.6;% ini
 rates = 0.010:0.001:0.2;% rec
%  rates = 0.010:0.001:0.6;% rel
dlogMBSdlogdG = zeros(size(rates,2),size(dG,2));
dlogBFdlogdG = zeros(size(rates,2),size(dG,2));


for r_idx = 1:1:size(rates,2)
    for idx = 1:1:size(dG,2)
        input_options.simulated_on      = true;
        input_options.EP_flag           = true;
        input_options.compute_TR        = false;
        input_options.attraction_coef   = 0.1; % K, E-P communication coefficient
        input_options.enhancer_index    = 25; % enhancer index number
        input_options.promoter_index    = 75; % promoter index number
        
        % Downstream: Transcriptional bursting
        input_options.k_on1_max         = 0.020; % burst initiation rate
        input_options.k_on1             = 0.002; % burst initiation rate
        input_options.k_recruitment_max = rates(r_idx); % burst polymerase recruitment rate
        input_options.k_recruitment     = 0.009; % burst polymerase recruitment rate
        input_options.k_release_max     = 0.025; % burst polymerase pause release rate
        input_options.k_release         = 0.008; % burst polymerase pause release rate
        input_options.k_off2            = 0.006;% off1 to off2
        input_options.k_on2             = 0.002;% off2 to off1
        input_options.k_off1            = 0.009;% rec to off1
        input_options.k_off3            = 0.002;% rec to off2
        params = ParametersBurst(input_options);
        params.EPGenomicDist = dG(idx);
        Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
            (abs(params.EPGenomicDist))+params.attraction_coef)^(-1));
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
        minlambda = min([params.k_on1,params.k_recruitment,params.k_release]);
        maxvelocity = (params.spring_const(1,1)/(abs(params.EPGenomicDist))...
            +params.attraction_coef)*d_EP(find(cumsum(P)>...
            99,1, 'first'))/params.friction_coef;
        mu = minlambda/maxvelocity;
       %%
        dFddG = 1./sqrt(2*pi).*(d_05./sqrt(params.diffusion_const.*params.friction_coef)).^3.*...
            sqrt((params.spring_const(1,1)./(params.EPGenomicDist)...
            +params.attraction_coef)).*exp(-0.5.*(d_05./sqrt(params.diffusion_const.*...
            params.friction_coef)).^2.*((params.spring_const(1,1)./(params.EPGenomicDist)...
            +params.attraction_coef))).*(-params.spring_const(1,1)./params.EPGenomicDist.^2);
        
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
        BS_meanV1 = (recruitment.*release./(recruitment.*off3+(release+off3).*off1));
        BS_meanV2 = sum((lambda_recruitment.*lambda_release./...
            (lambda_recruitment.*lambda_off3+...
            (lambda_release+lambda_off3).*lambda_off1)).*P.*dt);
        TON_meanV1 = (recruitment+release+off3)/(recruitment.*off3+(release+off3).*off1);
        TOFF_meanV1 = (off2+on2)./(on2.*on1);
        TON_meanV2 = sum((lambda_recruitment+lambda_release+lambda_off3)./...
            (lambda_recruitment.*lambda_off3+(lambda_release+lambda_off3).*lambda_off1).*P.*dt);
        TOFF_meanV2 = sum((lambda_off2+lambda_on2)./(lambda_on2.*lambda_on1).*P.*dt);
        CT_meanV1 = TON_meanV1 + TOFF_meanV1;
        CT_meanV2 = TON_meanV2 + TOFF_meanV2;
        
        dBSdlambdarel = (recruitment*off3*(recruitment+off1)/(recruitment*off3 + ...
            (release+off3)*off1)^2)*(params.k_release_max -params.k_release);
        dBSdlambdarec = (release*off1*(release+off3)/(recruitment*off3 + ...
            (release+off3)*off1)^2)*(params.k_recruitment_max -params.k_recruitment);
        dBSdlambda = dBSdlambdarel + dBSdlambdarec;
        UpperBS = (1./(1+mu)).*dBSdlambda+ (mu./(1+mu)).*(BS_max - BS_min);
        DownerBS = (1./(1+mu)).*BS_meanV1+ (mu./(1+mu)).*BS_meanV2;
        dlogMBSdlogdG(r_idx,idx) = params.EPGenomicDist.*dFddG.*(UpperBS./DownerBS);
        
        dCTdlambdaini = -(1+off2/on2)/(on1^2)*(params.k_on1_max -params.k_on1);
        dCTdlambdarec = (off1-off3)*(recruitment+off3)*(params.k_recruitment_max ...
            -params.k_recruitment)/(recruitment*off3 + (release+off3)*off1)^2;
        dCTdlambdarel = recruitment*(off3-off1)*(params.k_release_max...
            -params.k_release)/(recruitment*off3 + (release+off3)*off1)^2;
        dCTdlambda = dCTdlambdaini + dCTdlambdarec + dCTdlambdarel;            
        UpperCT = (1./(1+mu)).*dCTdlambda+ (mu./(1+mu)).*(CT_max - CT_min);
        DownerCT = (1./(1+mu)).*CT_meanV1+ (mu./(1+mu)).*CT_meanV2;
        dlogMCTdlogdG = params.EPGenomicDist.*dFddG.*(UpperCT./DownerCT);
        dlogBFdlogdG(r_idx,idx) = -dlogMCTdlogdG;
        
        
    end
    
    
end

%%

figure
S1 = surf(log10(dG),rates,dlogMBSdlogdG);
S1.EdgeColor = 'none';
hold on
S2 = surf(log10(dG),rates,dlogBFdlogdG,'FaceAlpha',0.8);
S2.EdgeColor = 'none';
box on



