% calculate mutal information Figure 5 & Figure S18
K = linspace(log10(0.01),log10(1),50);
k = linspace(log10(0.0001),log10(0.01),50);
IBSDS = zeros(size(K,2),size(k,2));
ICTDS = zeros(size(K,2),size(k,2));
IOFFDS= zeros(size(K,2),size(k,2));
IONDS = zeros(size(K,2),size(k,2));

for jdx = 1:1:50
    for idx = 1:1:size(K,2)
        input_options.EP_flag           = true;
        input_options.attraction_coef   = 10^(K(idx)); % K, E-P communication coefficient
        input_options.enhancer_index    = 25; % enhancer index number
        input_options.promoter_index    = 75; % promoter index number
        
        % Downstream: Transcriptional bursting
        input_options.k_on1_max         = 0.020; % burst initiation rate
        input_options.k_on1             = 10^(k(jdx));%0.002; % burst initiation rate
        input_options.k_recruitment_max = 0.024; % burst polymerase recruitment rate
        input_options.k_recruitment     = 0.009; % burst polymerase recruitment rate
        input_options.k_release_max     = 0.025; % burst polymerase pause release rate
        input_options.k_release         = 0.008; % burst polymerase pause release rate
        input_options.k_off2            = 0.006; % off1 to off2
        input_options.k_on2             = 0.002; % off2 to off1
        input_options.k_off1            = 0.009; % rec to off1
        input_options.k_off3            = 0.002; % rec to off2
        params = ParametersBurst(input_options);
        Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
            (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
        dt = 0.01;
        d_EP = 0.01:dt:3;
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
        maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
            params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
            99,1, 'first'))/params.friction_coef;
        mu = minlambda/maxvelocity;
        %% BSPDF
        Mmax = 50;
        m = 0:1:Mmax;
        theta =  (off3.*recruitment + (release+off3).*off1)./((recruitment+off1).*(release+off3));
        BSPDFV1 = theta.*(1-theta).^(m);
        
        m = 0:1:Mmax;
        theta = (lambda_off3.*lambda_recruitment + (lambda_release+lambda_off3).*...
            lambda_off1)./((lambda_recruitment+lambda_off1).*(lambda_release+lambda_off3));
        BS_prob = zeros(size(m,2),size(d_EP,2));
        for j = 1:size(m,2)
            BS_prob(j,:) = theta.*(1-theta).^(j-1);
        end
        BSPDFV2 = (BS_prob*(P.*dt)')';
        %% BS DS
        BSPmPdS = 1./(1+mu).*BSPDFV1'.*P +mu./(1+mu).*BS_prob.*P;
        BSPDF = 1./(1+mu).*BSPDFV1+mu./(1+mu).*BSPDFV2;
        IBSDS(idx,jdx) = sum(sum(BSPmPdS.*log2(BSPmPdS./(BSPDF'.*P)).*dt));
        
        %% CTPDF
        t = 0.01:0.01:10000;
        % ONPDFV1
        x1 = 0.5.*(-(recruitment+release+off1+off3)+sqrt((recruitment-release+off1-off3).^2 ...
            +4.*recruitment.*release));
        x2 = 0.5.*(-(recruitment+release+off1+off3)-sqrt((recruitment-release+off1-off3).^2 ...
            +4.*recruitment.*release));
        A1 = (recruitment+release+off3+x1).*x1;
        A2 = (recruitment+release+off3+x2).*x2;
        PONPDFV1 = -(A1.*exp(x1.*t)-A2.*exp(x2.*t))./(x1-x2);
        
        % OFFPDFV1
        y1 = 0.5.*(-(on2+on1+off2)+sqrt((off2+on1-on2).^2+4.*off2.*on2));
        y2 = 0.5.*(-(on2+on1+off2)-sqrt((off2+on1-on2).^2+4.*off2.*on2));
        B1 = (off2+on2+y1).*y1;
        B2 = (off2+on2+y2).*y2;
        POFFPDFV1 = -(B1.*exp(y1.*t)-B2.*exp(y2.*t))./(y1-y2);
        
        % CTPDFV1
        PCTPDFV1 = A1.*B1.*(exp(x1.*t)-exp(y1.*t))./((x1-y1).*(y1-y2).*(x1-x2))...
            +A2.*B2.*(exp(x2.*t)-exp(y2.*t))./((x2-y2).*(y1-y2).*(x1-x2))...
            -A2.*B1.*(exp(x2.*t)-exp(y1.*t))./((x2-y1).*(y1-y2).*(x1-x2))...
            -A1.*B2.*(exp(x1.*t)-exp(y2.*t))./((x1-y2).*(y1-y2).*(x1-x2));
        
        % ONPDFV2
        x1 = 0.5.*(-(lambda_recruitment+lambda_release+lambda_off1+lambda_off3)...
            +sqrt((lambda_recruitment-lambda_release+lambda_off1-lambda_off3).^2 ...
            +4.*lambda_recruitment.*lambda_release));
        x2 =  0.5.*(-(lambda_recruitment+lambda_release+lambda_off1+lambda_off3)...
            - sqrt((lambda_recruitment-lambda_release+lambda_off1-lambda_off3).^2 ...
            +4.*lambda_recruitment.*lambda_release));
        U_ON = (lambda_recruitment+lambda_release+lambda_off3+x1).*x1;
        V_ON = (lambda_recruitment+lambda_release+lambda_off3+x2).*x2;
        
        ON_prob = zeros(size(t,2),size(d_EP,2));
        for j = 1:size(t,2)
            ON_prob(j,:) = -(U_ON.*exp(x1.*t(j))-V_ON.*exp(x2.*t(j)))./(x1-x2);
        end
        PONPDFV2 = (ON_prob*(P.*dt)')';
        
        % OFFPDFV2
        y1 = 0.5.*(-(lambda_on2+lambda_on1+lambda_off2)+...
            sqrt((lambda_off2+lambda_on1-lambda_on2).^2+4.*lambda_on2.*lambda_off2));
        y2 = 0.5.*(-(lambda_on2+lambda_on1+lambda_off2)-...
            sqrt((lambda_off2+lambda_on1-lambda_on2).^2+4.*lambda_on2.*lambda_off2));
        U_OFF = (lambda_on2+lambda_off2+y1).*y1;
        V_OFF = (lambda_on2+lambda_off2+y2).*y2;
        OFF_prob = zeros(size(t,2),size(d_EP,2));
        for j = 1:size(t,2)
            OFF_prob(j,:) = -(U_OFF.*exp(y1.*t(j))-V_OFF.*exp(y2.*t(j)))./(y1-y2);
        end
        POFFPDFV2 = (OFF_prob*(P.*dt)')';
        
        % CTPDFV2
        CT_prob = zeros(size(t,2),size(d_EP,2));
        for j = 1:size(t,2)
            CT_prob(j,:) = U_ON.*U_OFF.*(exp(x1.*t(j))-exp(y1.*t(j)))./((x1-y1).*(y1-y2).*(x1-x2))...
                +V_ON.*V_OFF.*(exp(x2.*t(j))-exp(y2.*t(j)))./((x2-y2).*(y1-y2).*(x1-x2))...
                -V_ON.*U_OFF.*(exp(x2.*t(j))-exp(y1.*t(j)))./((x2-y1).*(y1-y2).*(x1-x2))...
                -U_ON.*V_OFF.*(exp(x1.*t(j))-exp(y2.*t(j)))./((x1-y2).*(y1-y2).*(x1-x2));
        end
        PCTPDFV2 = (CT_prob*(P.*dt)')';
        
        %% CT DS
        PONPdS = 1./(1+mu).*PONPDFV1'.*P +mu./(1+mu).*ON_prob.*P;
        PONPDF = 1./(1+mu).*PONPDFV1+mu./(1+mu).*PONPDFV2;
        IONDS(idx,jdx) = sum(sum(PONPdS.*log2(PONPdS./(PONPDF'.*P)).*dt.*0.01));
        
        POFFPdS = 1./(1+mu).*POFFPDFV1'.*P +mu./(1+mu).*OFF_prob.*P;
        POFFPDF = 1./(1+mu).*POFFPDFV1+mu./(1+mu).*POFFPDFV2;
        IOFFDS(idx,jdx) = sum(sum(POFFPdS.*log2(POFFPdS./(POFFPDF'.*P)).*dt.*0.01));
        
        PCTPdS = 1./(1+mu).*PCTPDFV1'.*P +mu./(1+mu).*CT_prob.*P;
        PCTPDF = 1./(1+mu).*PCTPDFV1+mu./(1+mu).*PCTPDFV2;
        ICTDS(idx,jdx) = sum(sum(PCTPdS.*log2(PCTPdS./(PCTPDF'.*P)).*dt.*0.01));
    end
end




%%
figure3 = figure;
set(figure3,'position',[300 400 280 190]);
imagesc(10.^k,10.^K,ICTDS);
axis xy
box on
set(gca,'TickLength',[0.02,0.025]);
axis square

figure4 = figure;
set(figure4,'position',[300 400 280 190]);
imagesc(10.^k,10.^K,IBSDS);
axis xy
box on
set(gca,'TickLength',[0.02,0.025]);
axis square


