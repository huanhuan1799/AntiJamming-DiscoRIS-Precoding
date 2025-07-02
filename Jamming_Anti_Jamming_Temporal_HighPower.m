%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main file for testing DIRS-based T-FPJ and the effect of anti-jamming precoding against T-FPJ
% The transmit power per LU is -2 dBm (high transmit power)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  clc

rng('shuffle'); % initiate the random number generators with a random seed
%% -------------------------------------------------------Parameter setting--------------------------------------------------------------------
Nr = 1; % number of LU's antennas
K = 12; % number of users
N = 16; % number of transmit antennas
Niter = 2000; % number of realizations in the Monte Carlo simulations
WW = 64;
HH = 32;
NN = WW*HH; % number of IRS elements

%Range of SNR values
PdBm = -2; % transmit power per LU (dBm)
P = (10.^(PdBm/10)); % linear scale (mW)
AJ_dB = -4; % jamming power of AJ (dB)
%%%%%%%%%%%%%%%%%%%%%%For Estimation
itemin = 1;
Tc = 6; % channel changing times at DT phase
itemax = Tc;
ite_Num = [itemin:1:itemax]; % number of feedback times for parctical estimation of ACA channel

% initialization
sumRateZF_wo = zeros(1,Niter);
sumRateZF_AJ = zeros(1,Niter);
 
sumRateZF_GC_Anti_P12_Frame = zeros(length(ite_Num),Niter);
sumRateZF_GC_Anti_P14_Frame = zeros(length(ite_Num),Niter);

sumRateZF_GC_Anti_P12 = zeros(1,Niter);
sumRateZF_GC_Anti_P14 = zeros(1,Niter);

sumRateZF_GC_TFPJ_P12 = zeros(1,Niter);
sumRateZF_GC_TFPJ_P14 = zeros(1,Niter);

%% Generate Random Reflecting Vector
clear Omg
b = 1;
Amp = [0.8, 1];
Omg = [pi/9,7*pi/6];
aaa = zeros(1,Niter);
bbb = zeros(1,Niter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
for m = 1:Niter
        m
%% Channel;Only DIRS, not random wireless channels
         eb = 10; eb22 = 1/(1+eb); eb11 = eb/(1+eb); eb1 = sqrt(eb11); eb2 = sqrt(eb22);
         [ Hd,Gr,Hr,lg,lr,ld,lAJ,noise] = J_generate_channel_nearfield(WW,HH,N,K,eb1,eb2);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% W/O RISs
W = (Hd'/(Hd*Hd')); % ZF precoding
for k = 1:K
    User_set = [1:K];  
    User_set(find(User_set == k)) = [];
    HD_age = zeros(N,N);
    for u = 1:length(User_set)
        HD_age = HD_age + Hd(User_set(u),:)'*Hd(User_set(u),:) ;
    end
    User_Sig = P *abs( (W(:,k)'./norm(W(:,k)','fro'))*...
        ( Hd(k,:)'*Hd(k,:)   ) *(W(:,k)./norm(W(:,k),'fro')) ); % effective signal power of LU k
    
    InterUser_Inter =  P *abs( (W(:,k)'./norm(W(:,k)','fro'))*...
        ( HD_age +  ( 1/P  ) *eye(N) )...
        *(W(:,k)./norm(W(:,k),'fro')) );                        % power of inter-user interference + AWGN
    sumRateZF_wo(m) = sumRateZF_wo(m) + log2( 1+User_Sig/(InterUser_Inter) ); % sum rate w/o RIS
end
           
clear wRPTZF_wo W User_Sig HD_age
 
%% Existing Active Jammer  
W = (Hd'/(Hd*Hd')); % ZF precoding
for k = 1:K
    User_set = [1:K]; InterUser_Inter = 0;
    User_set(find(User_set == k)) = [];
    for u = 1:length(User_set)
        InterUser_Inter = InterUser_Inter + (P)*abs(Hd(User_set(u),:)*(W(:,k)/norm(W(:,k),'fro')) )^2; % power of inter-user interference
    end
    User_Sig = (P )*abs(Hd(k,:)*(W(:,k)/norm(W(:,k),'fro')) )^2; % effective signal power of LU k
    sumRateZF_AJ(m) = sumRateZF_AJ(m) + log2( 1+User_Sig/(1+InterUser_Inter + lAJ(k)*10.^(-noise/10)*10^(AJ_dB/10) ) ); % sum rate when existing active jammer
end
                 
%% Generation Pha during the RPT and DT phases for T-FPJ
clear W User_Sig InterUser_Inter
% --------------------- RPT phase (case 1) --------------------- % 
HRPT_P14 =  Hd; % different from PFPJ, IRS is silent at the RPT phase of TFPJ
W_P14 = (HRPT_P14'/(HRPT_P14*HRPT_P14')) ; % ZF precoding
% -------------------------------------------------------------- %

% --------------------- RPT phase (case 2) --------------------- %
HRPT_P12 =  Hd;
W_P12 = (HRPT_P12'/(HRPT_P12*HRPT_P12')) ; % ZF precoding
% -------------------------------------------------------------- %

% -------------------------- DT phase -------------------------- %
theta_inbarDT_P14 = zeros(Tc,NN);
sumRateZF_Anti_P14_1 = zeros(1,Tc);
sumRateZF_TFPJ_P14_1 = zeros(1,Tc);

theta_inbarDT_P12 = zeros(Tc,NN);
sumRateZF_Anti_P12_1 = zeros(1,Tc);
sumRateZF_TFPJ_P12_1 = zeros(1,Tc);
for ic = 1:Tc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    betaDT_P14_temp = zeros(1,NN);
    phaDIRSDT_P14_Temp = zeros(1,NN);
    
    betaDT_P12_temp = zeros(1,NN);
    phaDIRSDT_P12_Temp = zeros(1,NN);
    for ri = 1:length(phaDIRSDT_P14_Temp)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iddele = randi(2^b,1,1)*randi(2^b,1,1) ;
        if iddele == 1 % 1/4 probability (case 1)
            phaDIRSDT_P14_Temp(ri) = Omg(1);
            betaDT_P14_temp(ri) = Amp(1); % phase-dependent amplitude
        else % 3/4 probability (case 1)
            phaDIRSDT_P14_Temp(ri) = Omg(2);
            betaDT_P14_temp(ri) = Amp(2); % phase-dependent amplitude
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iddele1 = randi(2^b,1,1);
        if iddele1 == 1 % 1/2 probability (case 2)
            phaDIRSDT_P12_Temp(ri) = Omg(1);
            betaDT_P12_temp(ri) = Amp(1); % phase-dependent amplitude
        else % 1/2 probability (case 2)
            phaDIRSDT_P12_Temp(ri) = Omg(2);
            betaDT_P12_temp(ri) = Amp(2); % phase-dependent amplitude
        end
    end
    theta_inbarDT_P14(ic,:) = betaDT_P14_temp.*exp(1j.*(phaDIRSDT_P14_Temp)); % random reflecting vector (time-changing) at DT phase (jamming), case 1
    theta_inbarDT_P12(ic,:) = betaDT_P12_temp.*exp(1j.*(phaDIRSDT_P12_Temp)); % random reflecting vector (time-changing) at DT phase (jamming), case 2
% -------------------------------------------------------------- %


%% Anti-Jamming Precoding
clear HD_DT_EachTime_P14 HD_DT_EachTime_P12
GC_Stacha_Acc_P12 = lr*lg*NN* 0.82 *(10.^(-noise/10)); % α^{_} = 0.91 case 1 
GC_Stacha_Acc_P14 = lr*lg*NN* 0.91 *(10.^(-noise/10)); % α^{_} = 0.82 case 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HD_DT_EachTime_P14 = Hr*diag( theta_inbarDT_P14(ic,:) )*Gr + Hd; % equivalent channel at different time, case 1 
HD_DT_EachTime_P12 = Hr*diag( theta_inbarDT_P12(ic,:) )*Gr + Hd; % equivalent channel at different time, case 2

sumRateZF_Anti_P14_11 = 0;
sumRateZF_Anti_P12_11 = 0;

W_Anti_P14 = zeros(N,K);
W_Anti_P12 = zeros(N,K);

% --------------------------- case 1 --------------------------- %
for k = 1:K
    clear Ak_The_P14
    User_set = [1:K];
    User_set(find(User_set == k)) = [];
    HD_age_The_P14 = zeros(N,N);
    HD_age_P14 = zeros(N,N);

    for u = 1:length(User_set)
        HD_age_The_P14  =   HD_age_The_P14 + HRPT_P14(User_set(u),:)'*HRPT_P14(User_set(u),:)+GC_Stacha_Acc_P14(User_set(u))*eye(N);
        HD_age_P14   =  HD_age_P14 + HD_DT_EachTime_P14( User_set(u),:)'*HD_DT_EachTime_P14( User_set(u),:);
    end
    Ak_The_P14 = inv( HD_age_The_P14 + (1/P) *eye(N) ) * ...
        ( (HRPT_P14(k,:))'*HRPT_P14(k,:) +  GC_Stacha_Acc_P14(k)*eye(N) ); % A_k case 1

    [VV,DD] = eig(Ak_The_P14);
    ak_the_P14 = find( abs(diag(DD)) == max(abs(diag(DD))) );
    W_Anti_P14(:,k) = VV(:,ak_the_P14); % find the eigenvector related to the max singular value

    User_Sig_The = P *abs( ...
        (W_Anti_P14(:,k)'./norm(W_Anti_P14(:,k)','fro'))*( (HD_DT_EachTime_P14(k,:))'*HD_DT_EachTime_P14(k,:) )...
        *(W_Anti_P14(:,k)./norm(W_Anti_P14(:,k),'fro')) ); % effective signal power

    InterUser_Inter_The =  P*abs( (W_Anti_P14(:,k)'./norm(W_Anti_P14(:,k)','fro'))*...
        ( HD_age_P14 +  ( 1/P ) *eye(N) )...
        *(W_Anti_P14(:,k)./norm(W_Anti_P14(:,k),'fro')) ); % power of inter-user interference and AWGN

    sumRateZF_Anti_P14_11  = sumRateZF_Anti_P14_11  + log2( 1+User_Sig_The/(InterUser_Inter_The) ); % sum rate when anti-jamming precoding is utilized
    clear  W_Anti  User_Sig_The InterUser_Inter_The HD_age_P14 Ak_The_P14 Ak_The_P14 ak_the_P14
end
sumRateZF_Anti_P14_1(ic) = sumRateZF_Anti_P14_11;
% -------------------------------------------------------------- %

% --------------------------- case 2 --------------------------- %
for k = 1:K
    clear Ak_The_P12
    User_set = [1:K];
    User_set(find(User_set == k)) = [];
    HD_age_The_P12 = zeros(N,N);
    HD_age_P12 = zeros(N,N);

    for u = 1:length(User_set)
        HD_age_The_P12  =   HD_age_The_P12 + HRPT_P12(User_set(u),:)'*...
            HRPT_P12(User_set(u),:)+GC_Stacha_Acc_P12(User_set(u))*eye(N);

        HD_age_P12   =  HD_age_P12 + HD_DT_EachTime_P12( User_set(u),:)'*HD_DT_EachTime_P12( User_set(u),:);
    end
    Ak_The_P12 = inv( HD_age_The_P12 + (1/P) *eye(N) ) * ...
        ( (HRPT_P12(k,:))'*HRPT_P12(k,:) +  GC_Stacha_Acc_P12(k)*eye(N) ); % A_k case 2

    [VV,DD] = eig(Ak_The_P12);
    ak_the_P12 = find( abs(diag(DD)) == max( abs(diag(DD)) ) );
    W_Anti_P12(:,k) = VV(:,ak_the_P12);

    User_Sig_The = P *abs( ...
        (W_Anti_P12(:,k)'./norm(W_Anti_P12(:,k)','fro'))*( (HD_DT_EachTime_P12(k,:))'*HD_DT_EachTime_P12(k,:) )...
        *(W_Anti_P12(:,k)./norm(W_Anti_P12(:,k),'fro')) ); % effective signal power

    InterUser_Inter_The =  P*abs( (W_Anti_P12(:,k)'./norm(W_Anti_P12(:,k)','fro'))*...
        ( HD_age_P12 +  ( 1/P ) *eye(N) )...
        *(W_Anti_P12(:,k)./norm(W_Anti_P12(:,k),'fro')) ); % power of inter-user interference and AWGN

    sumRateZF_Anti_P12_11  = sumRateZF_Anti_P12_11  + log2( 1+User_Sig_The/(InterUser_Inter_The) ); % sum rate when anti-jamming precoding is utilized
    clear  W_Anti  User_Sig_The InterUser_Inter_The HD_age_P12 Ak_The_P12 Ak_The_P12 ak_the_P12
end
sumRateZF_Anti_P12_1(ic) = sumRateZF_Anti_P12_11;   
% -------------------------------------------------------------- %   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Temporal Fully Passive Jammer, case 1: 1/4
sumRateZF_TFPJ_P14_11 = 0;  
    clear HD_DT_EachTime_P14
    HD_DT_EachTime_P14 = Hr*diag( theta_inbarDT_P14(ic,:) )*Gr + Hd;
for k = 1:K
    User_set = [1:K];
    User_set(find(User_set == k)) = [];
    HD_age_P14 = zeros(N,N);
    for u = 1:length(User_set)
        HD_age_P14   =  HD_age_P14 + HD_DT_EachTime_P14( User_set(u),:)'*HD_DT_EachTime_P14( User_set(u),:);
    end
    User_Sig = P *abs( (W_P14(:,k)'./norm(W_P14(:,k)','fro'))*...
        ( HD_DT_EachTime_P14(k,:)'*HD_DT_EachTime_P14(k,:) )...
        *(W_P14(:,k)./norm(W_P14(:,k),'fro')) ); % effective signal power
     
    InterUser_Inter =  P *abs( (W_P14(:,k)'./norm(W_P14(:,k)','fro'))*...
        ( HD_age_P14 +  ( 1/P ) *eye(N) ) ...
        *(W_P14(:,k)./norm(W_P14(:,k),'fro')) ) ; % power of inter-user interference (including ACA interference) and AWGN
     
    sumRateZF_TFPJ_P14_11 = sumRateZF_TFPJ_P14_11 + log2( 1+User_Sig/(InterUser_Inter) ); % sum rate of TFPJ
    clear InterUser_Inter User_Sig
end  
sumRateZF_TFPJ_P14_1(ic) = sumRateZF_TFPJ_P14_11;

%% Temporal Fully Passive Jammer, case 2: 1/2
sumRateZF_TFPJ_P12_11 = 0;
clear HD_DT_EachTime_P12
HD_DT_EachTime_P12 = Hr*diag( theta_inbarDT_P12(ic,:) )*Gr + Hd;
for k = 1:K
    User_set = [1:K];
    User_set(find(User_set == k)) = [];
    HD_age_P12 = zeros(N,N);
    for u = 1:length(User_set)
        HD_age_P12   =  HD_age_P12 + HD_DT_EachTime_P12( User_set(u),:)'*HD_DT_EachTime_P12( User_set(u),:);
    end
    User_Sig = P *abs( (W_P12(:,k)'./norm(W_P12(:,k)','fro'))*...
        ( HD_DT_EachTime_P12(k,:)'*HD_DT_EachTime_P12(k,:) )...
        *(W_P12(:,k)./norm(W_P12(:,k),'fro')) ); % effective signal power

    InterUser_Inter =  P *abs( (W_P12(:,k)'./norm(W_P12(:,k)','fro'))*...
        ( HD_age_P12 +  ( 1/P ) *eye(N) ) ...
        *(W_P12(:,k)./norm(W_P12(:,k),'fro')) ); % power of inter-user interference (including ACA interference) and AWGN

    sumRateZF_TFPJ_P12_11 = sumRateZF_TFPJ_P12_11 + log2( 1+User_Sig/(InterUser_Inter) ); % sum rate of TFPJ
    clear InterUser_Inter User_Sig
end  
sumRateZF_TFPJ_P12_1(ic) = sumRateZF_TFPJ_P12_11;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end
sumRateZF_GC_Anti_P14(m) = sum(sumRateZF_Anti_P14_1)/Tc;
sumRateZF_GC_Anti_P12(m) = sum(sumRateZF_Anti_P12_1)/Tc;

sumRateZF_GC_TFPJ_P14(m) = sum(sumRateZF_TFPJ_P14_1)/Tc;
sumRateZF_GC_TFPJ_P12(m) = sum(sumRateZF_TFPJ_P12_1)/Tc;
 
%% Practical estimation of Statistical Characteristics of ACA channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it = 1:length(ite_Num)
    CurrMax_Ite = ite_Num(it); % Max Number of feedback times 
    Est_s_P14 = zeros(K,CurrMax_Ite);
    %     H_DT_all_P14 = zeros(K,N,CurrMax_Ite);
    sumRateZF_Anti_P14_temp = 0;
    sumRateZF_Anti_P14_temp1 = zeros(1,Tc);
    
    Est_s_P12 = zeros(K,CurrMax_Ite);
    %     H_DT_all_P12 = zeros(K,N,CurrMax_Ite);
    sumRateZF_Anti_P12_temp = 0;
    sumRateZF_Anti_P12_temp1 = zeros(1,Tc);
    for is = 1:CurrMax_Ite
        %% Estimate Sta Vau
        HD_DT_ACA_P14_power = zeros(1,K);
        HD_DT_ACA_P12_power = zeros(1,K);
        for k = 1:K
            for iis = 1:is
                clear H_DT_P14_temp H_DT_P12_temp
                %                 H_DT_all_P14(:,:,is) = Hr*diag(theta_inbarDT_P14(iis,:))*Gr + Hd;
                %                 H_DT_all_P12(:,:,is) = Hr*diag(theta_inbarDT_P12(iis,:))*Gr + Hd;
                H_DT_P14_temp(:,:) = Hr*diag(theta_inbarDT_P14(iis,:))*Gr + Hd;
                H_DT_P12_temp(:,:) = Hr*diag(theta_inbarDT_P12(iis,:))*Gr + Hd;
                HD_DT_ACA_P14_power(k) = HD_DT_ACA_P14_power(k) + norm(H_DT_P14_temp(k,:),'fro')^2; % feedback received power values
                HD_DT_ACA_P12_power(k) = HD_DT_ACA_P12_power(k) + norm(H_DT_P12_temp(k,:),'fro')^2;
            end
        end
        HD_DT_ACA_P14_power = HD_DT_ACA_P14_power./is;
        HD_DT_ACA_P12_power = HD_DT_ACA_P12_power./is;
        for k = 1:K
            Est_s_P14(k,is) =  abs(  HD_DT_ACA_P14_power(k) - norm(HRPT_P14(k,:),'fro')^2 )/N; % the s-th estimated statistical characteristic
            Est_s_P12(k,is) =  abs(  HD_DT_ACA_P12_power(k) - norm(HRPT_P12(k,:),'fro')^2 )/N;
        end

        %% Proposed Anti-Jamming in This Trans w/ the sth est vua
        W_Anti_P14 = zeros(N,K);
        W_Anti_P12 = zeros(N,K);
        clear HD_DT_EachTime_P14 HD_DT_EachTime_P12
        %HD_Aging_EachTime = Hr*diag( theta_inbarDT_P14(ic,:)-theta_inbarRPT_P14 )*Gr + Hd;
        HD_DT_EachTime_P14 = Hr*diag( theta_inbarDT_P14(is,:) )*Gr + Hd;
        HD_DT_EachTime_P12 = Hr*diag( theta_inbarDT_P12(is,:) )*Gr + Hd;
        for k = 1:K
            clear Ak_The Ak_Simu
            User_set = [1:K];
            User_set(find(User_set == k)) = [];
            HD_age_The_P14 = zeros(N,N);
            HD_age_P14 = zeros(N,N);

            for u = 1:length(User_set)
                HD_age_The_P14  =   HD_age_The_P14 + HRPT_P14(User_set(u),:)'*HRPT_P14(User_set(u),:)+Est_s_P14(User_set(u),is)*eye(N); % estimated
                HD_age_P14   =  HD_age_P14 + HD_DT_EachTime_P14( User_set(u),:)'*HD_DT_EachTime_P14( User_set(u),:); % actual
            end
            Ak_The_P14 = inv( HD_age_The_P14 + (1/P) *eye(N) ) * ( (HRPT_P14(k,:))'*HRPT_P14(k,:) +  Est_s_P14(k,is)*eye(N) );

            [VV,DD] = eig(Ak_The_P14);
            ak_the_P14 = find( diag(DD) == max(diag(DD)) );
            W_Anti_P14(:,k) = VV(:,ak_the_P14);

            User_Sig_The = P *abs( ...
                (W_Anti_P14(:,k)'./norm(W_Anti_P14(:,k)','fro'))*( (HD_DT_EachTime_P14(k,:))'*HD_DT_EachTime_P14(k,:) )...
                *(W_Anti_P14(:,k)./norm(W_Anti_P14(:,k),'fro')) );

            InterUser_Inter_The =  P*abs( (W_Anti_P14(:,k)'./norm(W_Anti_P14(:,k)','fro'))*...
                ( HD_age_P14 +  ( 1/P ) *eye(N) )...
                *(W_Anti_P14(:,k)./norm(W_Anti_P14(:,k),'fro')) );

            sumRateZF_Anti_P14_temp1(is) = sumRateZF_Anti_P14_temp1(is) + log2( 1+User_Sig_The/(InterUser_Inter_The) ); % sum rate when anti jamming precoding is utlized based on estimated ACA channels
            clear W_Anti  User_Sig_The InterUser_Inter_The User_Sig_Simu InterUser_Inter_Simu  HD_age_P14   Ak_The_P14 Ak_The_P14 ak_the_P14
        end
        % -------------------------------------------------------------- %

        % --------------------------- case 2 --------------------------- %
        for k = 1:K
            clear Ak_The_P12
            User_set = [1:K];
            User_set(find(User_set == k)) = [];
            HD_age_The_P12 = zeros(N,N);
            HD_age_P12 = zeros(N,N);

            for u = 1:length(User_set)
                HD_age_The_P12  =   HD_age_The_P12 + HRPT_P12(User_set(u),:)'*HRPT_P12(User_set(u),:)+Est_s_P12(User_set(u),is)*eye(N); % estimated
                HD_age_P12   =  HD_age_P12 + HD_DT_EachTime_P12( User_set(u),:)'*HD_DT_EachTime_P12( User_set(u),:); % actual
            end
            Ak_The_P12 = inv( HD_age_The_P12 + (1/P) *eye(N) ) * ( (HRPT_P12(k,:))'*HRPT_P12(k,:) +  Est_s_P12(k,is)*eye(N) );

            [VV,DD] = eig(Ak_The_P12);
            ak_the_P12 = find( diag(DD) == max(diag(DD)) );
            W_Anti_P12(:,k) = VV(:,ak_the_P12);

            User_Sig_The = P *abs( ...
                (W_Anti_P12(:,k)'./norm(W_Anti_P12(:,k)','fro'))*( (HD_DT_EachTime_P12(k,:))'*HD_DT_EachTime_P12(k,:) )...
                *(W_Anti_P12(:,k)./norm(W_Anti_P12(:,k),'fro')) );

            InterUser_Inter_The =  P*abs( (W_Anti_P12(:,k)'./norm(W_Anti_P12(:,k)','fro'))*...
                ( HD_age_P12 +  ( 1/P ) *eye(N) )...
                *(W_Anti_P12(:,k)./norm(W_Anti_P12(:,k),'fro')) );

            sumRateZF_Anti_P12_temp1(is) = sumRateZF_Anti_P12_temp1(is) + log2( 1+User_Sig_The/(InterUser_Inter_The) ); % sum rate when anti jamming precoding is utlized based on estimated ACA channels
            clear W_Anti  User_Sig_The InterUser_Inter_The User_Sig_Simu InterUser_Inter_Simu  HD_age_P12  Ak_The_P12 Ak_The_P12 ak_the_P12
        end
        % -------------------------------------------------------------- %
    end

    %% the last Tc-b bits
    for itt = 1:Tc - CurrMax_Ite
        Ite = itt + CurrMax_Ite;
        clear HD_DT_EachTime_P14 HD_DT_EachTime_P12
        HD_DT_EachTime_P14 = Hr*diag( theta_inbarDT_P14(Ite,:) )*Gr + Hd;
        HD_DT_EachTime_P12 = Hr*diag( theta_inbarDT_P12(Ite,:) )*Gr + Hd;

        % --------------------------- case 1 --------------------------- %
        for k = 1:K
            clear Ak_The_P14
            User_set = [1:K];
            User_set(find(User_set == k)) = [];
            HD_age_The_P14 = zeros(N,N);
            HD_age_P14 = zeros(N,N);

            for u = 1:length(User_set)
                HD_age_The_P14  =   HD_age_The_P14 + HRPT_P14(User_set(u),:)'*HRPT_P14(User_set(u),:)+Est_s_P14(User_set(u),CurrMax_Ite)*eye(N); % estimated
                HD_age_P14   =  HD_age_P14 + HD_DT_EachTime_P14( User_set(u),:)'*HD_DT_EachTime_P14( User_set(u),:); % actual
            end
            Ak_The_P14 = inv( HD_age_The_P14 + (1/P) *eye(N) ) * ( (HRPT_P14(k,:))'*HRPT_P14(k,:) +  Est_s_P14(k,CurrMax_Ite)*eye(N) );

            [VV,DD] = eig(Ak_The_P14);
            ak_the_P14 = find( diag(DD) == max(diag(DD)) );
            W_Anti_P14(:,k) = VV(:,ak_the_P14);

            User_Sig_The = P *abs( ...
                (W_Anti_P14(:,k)'./norm(W_Anti_P14(:,k)','fro'))*( (HD_DT_EachTime_P14(k,:))'*HD_DT_EachTime_P14(k,:) )...
                *(W_Anti_P14(:,k)./norm(W_Anti_P14(:,k),'fro')) );

            InterUser_Inter_The =  P*abs( (W_Anti_P14(:,k)'./norm(W_Anti_P14(:,k)','fro'))*...
                ( HD_age_P14 +  ( 1/P ) *eye(N) )...
                *(W_Anti_P14(:,k)./norm(W_Anti_P14(:,k),'fro')) );

            sumRateZF_Anti_P14_temp1(Ite) = sumRateZF_Anti_P14_temp1(Ite) + log2( 1+User_Sig_The/(InterUser_Inter_The) ); % sum rate when anti jamming precoding is utlized based on estimated ACA channels
            clear W_Anti  User_Sig_The InterUser_Inter_The User_Sig_Simu InterUser_Inter_Simu  HD_age_P14 Ak_The_P14 Ak_The_P14 ak_the_P14
        end
        % -------------------------------------------------------------- %

        % --------------------------- case 2 --------------------------- %
        for k = 1:K
            clear Ak_The_P12
            User_set = [1:K];
            User_set(find(User_set == k)) = [];
            HD_age_The_P12 = zeros(N,N);
            HD_age_P12 = zeros(N,N);

            for u = 1:length(User_set)
                HD_age_The_P12  =   HD_age_The_P12 + HRPT_P12(User_set(u),:)'*HRPT_P12(User_set(u),:)+Est_s_P12(User_set(u),CurrMax_Ite)*eye(N);  % estimated
                HD_age_P12   =  HD_age_P12 + HD_DT_EachTime_P12( User_set(u),:)'*HD_DT_EachTime_P12( User_set(u),:); % actual
            end
            Ak_The_P12 = inv( HD_age_The_P12 + (1/P) *eye(N) ) * ( (HRPT_P12(k,:))'*HRPT_P12(k,:) +  Est_s_P12(k,CurrMax_Ite)*eye(N) );

            [VV,DD] = eig(Ak_The_P12);
            ak_the_P12 = find( diag(DD) == max(diag(DD)) );
            W_Anti_P12(:,k) = VV(:,ak_the_P12);

            User_Sig_The = P *abs( ...
                (W_Anti_P12(:,k)'./norm(W_Anti_P12(:,k)','fro'))*( (HD_DT_EachTime_P12(k,:))'*HD_DT_EachTime_P12(k,:) )...
                *(W_Anti_P12(:,k)./norm(W_Anti_P12(:,k),'fro')) );

            InterUser_Inter_The =  P*abs( (W_Anti_P12(:,k)'./norm(W_Anti_P12(:,k)','fro'))*...
                ( HD_age_P12 + ( 1/P ) *eye(N) )...
                *(W_Anti_P12(:,k)./norm(W_Anti_P12(:,k),'fro')) );

            sumRateZF_Anti_P12_temp1(Ite) = sumRateZF_Anti_P12_temp1(Ite) + log2( 1+User_Sig_The/(InterUser_Inter_The) );  % sum rate when anti jamming precoding is utlized based on estimated ACA channels
            clear W_Anti  User_Sig_The InterUser_Inter_The User_Sig_Simu InterUser_Inter_Simu  HD_age_P14 Ak_The_P14 Ak_The_P14 ak_the_P14
        end
        % -------------------------------------------------------------- %
    end

    sumRateZF_GC_Anti_P14_Frame(it,m) = sum(sumRateZF_Anti_P14_temp1)/Tc; % time-average
    sumRateZF_GC_Anti_P12_Frame(it,m) = sum(sumRateZF_Anti_P12_temp1)/Tc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
end
        

%Plot simulation results
close all
figure; hold on; box on;
%% MU
yyaxis left
plot(ite_Num,(sum(sumRateZF_wo)/K/Niter)*ones(1,length(ite_Num)) ,'k-*','LineWidth',1.8,'MarkerSize',8); hold on 

plot(ite_Num,(sum(sumRateZF_GC_Anti_P14)/K/Niter)*ones(1,length(ite_Num)) ,'r-.v','LineWidth',1.8,'MarkerSize',8); hold on 
plot(ite_Num,mean(sumRateZF_GC_Anti_P14_Frame(:,:),2)/K,'cp','LineWidth',1.8,'MarkerSize',8);hold on 

plot(ite_Num,(sum(sumRateZF_GC_Anti_P12)/K/Niter)*ones(1,length(ite_Num)) ,'b-.o','LineWidth',1.8,'MarkerSize',8);hold on 
plot(ite_Num,mean(sumRateZF_GC_Anti_P12_Frame(:,:),2)/K,'c<','LineWidth',1.8,'MarkerSize',8);hold on 

plot(ite_Num,(sum(sumRateZF_GC_TFPJ_P14)/K/Niter)*ones(1,length(ite_Num)) ,'r--s','LineWidth',1.8,'MarkerSize',8);hold on 
plot(ite_Num,(sum(sumRateZF_GC_TFPJ_P12)/K/Niter)*ones(1,length(ite_Num)) ,'b--d','LineWidth',1.8,'MarkerSize',8);hold on 

plot(ite_Num,(sum(sumRateZF_AJ)/K/Niter)*ones(1,length(ite_Num)) ,'g-h','LineWidth',1.8,'MarkerSize',8);
ylabel('Ergodic Rate Per LU [bits/symbol/user]');
% legend('W/O Jamming','Proposed AJP & C1','Proposed Est & C1','W/ TFPJ & C1',...
%     'Proposed AJP & C2','Proposed Est & C2','W/ TFPJ & C2','W/ {P_J} = -4 dBm','Location','Best');

yyaxis right
clear yy_temp yyy_temp yyyy_temp zz_temp zzz_temp zzzz_temp
yy_temp = mean(sumRateZF_GC_Anti_P14_Frame(:,:),2)/K;
yyy_temp = (sum(sumRateZF_GC_TFPJ_P14)/K/Niter)*ones(1,length(ite_Num));
yyyy_temp = ((yy_temp.' - yyy_temp)./yyy_temp)*100;

zz_temp = mean(sumRateZF_GC_Anti_P12_Frame(:,:),2)/K;
zzz_temp = (sum(sumRateZF_GC_TFPJ_P12)/K/Niter)*ones(1,length(ite_Num));
zzzz_temp = ((zz_temp.' - zzz_temp)./zzz_temp)*100;

plot(ite_Num ,yyyy_temp,'m:x','LineWidth',1.8,'MarkerSize',8);
hold on
plot(ite_Num ,zzzz_temp,'m:+','LineWidth',1.8,'MarkerSize',8);

ylabel('Anti-Jamming Gain Per LU [%]');
grid off
xlabel('Number of Feedback')
axis([itemin,itemax,-inf,inf])

legend('W/O Jamming','Proposed AJP & C1','Proposed Est & C1',...
    'Proposed AJP & C2','Proposed Est & C2','W/ T-FPJ & C1','W/ T-FPJ & C2','W/ {P_J} = -4 dBm',...
    'Gain w/ C1','Gain w/ C2','Location','Best');