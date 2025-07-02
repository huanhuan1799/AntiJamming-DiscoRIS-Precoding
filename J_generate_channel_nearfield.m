function [Hd,G,Hr,path_g,path_r,path_d,path_AJ,noise] = J_generate_channel_nearfield(WW,HH,M,K,eb1,eb2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for gennerating all channels

% WW: number of columns of DIRS
% HH: number of rows of DIRS
% M: number of transmit antennas of AP
% K: number of LUs
% eb1: LOS parameter
% eb2: NLOS parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting
N = WW*HH; % number of DIRS elements
AP = [0,0,5]; % location of AP
xxx= 2;
IRS = [xxx,0,5]; % location of DIRS
AJ = IRS; % location of AJ
noise = -170+10*log10(180*1e3); % AWGN power (dBm) with frequency band 180K Hz

%% location of LUs (random in a circle)
Rai = 20;
Tha = 2*pi*(rand(1,K));
rad(1:K) = Rai*rand(1,K);
Dar1  = 0;
Dar2 = 180;

for k = 1:K
    Tha(k) = 2*pi*rand(1);
    rad(k) = Rai*rand(1);
    Pt(k,:) = [Dar1+rad(k)*cos(Tha(k)),Dar2+rad(k)*sin(Tha(k)),0];
end

%% large-scale path loss of the AJ-LU link
for k = 1:K
    dAJ(k) = sqrt(sum(abs(Pt(k,:)-AJ).^2)); % distance between AJ and LU k
    LAJ(k) = 32.6 + 36.7*log10(dAJ(k)); % large-scale path loss (dB)
    path_AJ(k) = 10.^((-LAJ(k))/10); % large-scale path loss
    pAJ(k) = sqrt(path_AJ(k)*10.^(-noise/10)); % normalization
end

%% large-scale path loss of the AP-IRS link
d_g = sqrt(sum(abs(AP-IRS).^2));% distance between AP and IRS
Lg = 35.6 + 22*log10(d_g); % large-scale path loss (dB)
path_g = 10.^((-Lg)/10); % large-scale path loss
pg = sqrt(path_g*10.^(-noise/10)); % normalization

%% channel between AP and IRS (near filed)
G_LOS = zeros(N,M);
G_NLOS = sqrt(1/2).*(randn(N,M)+1j.*randn(N,M));
for  w = 1:WW
    for h = 1:HH
        for m = 1:M           
%             D_m_hw = norm(IRS_hw_ele - AP_m_ant,'fro'); % distance between the h-th row and w-th column IRS element and the m-th AP antenna
%             D_m = norm(IRS - AP_m_ant,'fro'); % distance between the centre of IRS and the m-th AP antenna
%             Ddiff = D_m_hw - D_m;          
            Ddiff = norm([(m-1)*0.05/2, 0, 2] - [xxx+(w-1)*0.05/2, (h-1)*0.05/2, 2],'fro'); % distance difference, Ddiff = D_m_hw - D_m;
            r = (w-1)*HH+h;
            G_LOS(r,m) = exp(-1j*(2*pi/0.05)*Ddiff);
        end
    end
end
G = eb1.*G_LOS+eb2.*G_NLOS;
G = pg.*G;

%% channel between IRS and LUs
Hr_sig = sqrt(1/2).*(randn(K,N)+1j.*randn(K,N));
for k = 1:K
    %% channel Hr_w
    dr(k) = sqrt(sum(abs(Pt(k,:)-IRS).^2)); % distance between IRS and LU k
    Lr(k) = 32.6 + 36.7*log10(dr(k)); % large-scale path loss (dB)
    path_r(k) = 10.^((-Lr(k))/10); % large-scale path loss
    pr(k) = sqrt(path_r(k));
    Hr(k,:) = pr(k).*Hr_sig(k,:);
end

%% direct channel between AP and LUs
for k = 1:K
    %% channel Hd_w
    dk(k) = sqrt(sum(abs(Pt(k,:)-AP).^2)); % distance between AP and LU k
    Ld(k) = 32.6 + 36.7*log10(dk(k)); % large-scale path loss (dB)
    path_d(k) = 10.^((-Ld(k))/10); % large-scale path loss
    pd(k) = sqrt(path_d(k)*10.^(-noise/10)); % normalization
    Hd(k,:) = sqrt(1/2).*(randn(1,M)+1j.*randn(1,M));
    Hd(k,:) = pd(k).*Hd(k,:);
end

