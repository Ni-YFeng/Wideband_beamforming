clear all;

%% 1. 系统参数初始化
% initialization
C= 4;%簇模型的簇数
R = 10;%每个簇有10条子路径
N_t = 256;  %antenna number
N_s = 4;%用户的天线数
N_r = 4;%用户的传输流数
N_RF = 4;%基站射频链路数

fc=10e10;%中心载频
fs=0.1*fc;%带宽
tmax=25e-9;%最大中心时延
sigma_t = 0.06e-9;%每个簇内的时延扩展
num=64;   %子载波数
SNR_dB= -20:5:30;%信噪比
SNR = 10.^(SNR_dB/10.);

N_iter = 50;%迭代次数
sigma_bs = 5 * pi / 180;%基站端角度扩展
sigma_mt = 5 * pi / 180;%用户端角度扩展

K = 16;%每个射频链路连接的延时器数量
P = N_t/K;%每个延时器连接的移相器数量 (P=256/16=16)

%% 2. 构造天线阵列的导引矢量字典 (用于OMP算法)
 theta_list = linspace(-1,1,512)*pi/2;% 将空间角度划分为512个网格
 At = zeros(length(theta_list), N_t);
 for j = 1:length(theta_list)
     % 根据公式(4)生成中心频点下的阵列响应向量
    At(j,:)=1/sqrt(N_t)*exp(-1j*[0:N_t-1]*pi*sin(theta_list(j)));
 end
 
 %% 3. 蒙特卡洛主循环

Rate=zeros(1,length(SNR));
Rate1=zeros(1,length(SNR));
Rate2=zeros(1,length(SNR));
Rate3=zeros(1,length(SNR));
Rate4=zeros(1,length(SNR));
Rate5=zeros(1,length(SNR));
Rate6=zeros(1,length(SNR));
Rate7=zeros(1,length(SNR));

temp0=0;temp1=0;temp2=0;temp3=0;temp4=0;temp5=0; temp6 = 0;temp7 = 0;

for i =1:N_iter
    tic()
    % 3.1 生成太赫兹宽带多径/簇信道
    % H为包含所有子载波的宽带信道，hc为中心频点的信道
     [H,hc] = WidebandChannel_Multi(N_t,N_r,C,R,fc,fs,tmax,num,sigma_bs, sigma_mt,sigma_t);

     % 3.2 计算最优无约束全数字预编码矩阵 (F_opt)
     F_opt = zeros(N_t, N_s, num);
     for j=1:num
         [~,~,V]=svd(H(:,:,j));% 对每个子载波的信道进行SVD分解
         F_opt(:,:,j) = V(:,1:N_s);% 取最大的N_s个右奇异向量作为最优预编码
     end
     [~,~,V]=svd(hc);
    V = V(:,1:N_s);
    F_optc = V;% 中心频点的最优预编码

%% 4. 各算法的速率计算与对比  
    %% 算法1：最优全数字预编码 (理论上限)
     F_max = zeros(num,N_t);
    
     for i_snr = 1:length(SNR)
     for j=1:num
         F = F_opt(:,:,j);
         HH = H(:,:,j)*F;
         % 香农容量公式：log2(det(I + H*F*F'*H'*SNR))
         Rate(i_snr) = Rate(i_snr) + real(log2(det(eye(N_s)+HH'*HH*SNR(i_snr))));
     end
     end
  
     
     %% 算法2：OMP (基于窄带的空域稀疏预编码, 文献[19])
     % 仅根据中心频点信道 hc 提取空间角度，生成频率无关的模拟波束
     [F1_temp,F_BB,th]=spatially_sparse_precoding(N_RF,N_s,F_optc,At.');
     theta = -theta_list(th);
     for i_snr = 1:length(SNR)
     for j=1:num
        Fbb1 = F1_temp\F_opt(:,:,j);
        Fbb1 = Fbb1*sqrt(N_s)/norm(F1_temp*Fbb1,'F');
        F1 = F1_temp*Fbb1;
        HH1 = H(:,:,j)*F1;
        Rate1(i_snr) = Rate1(i_snr) + real(log2(det(eye(N_s)+HH1'*HH1*SNR(i_snr))));
     end
     end
     
     
    %% 算法3：空时分组编码 (STBC, 文献[8])
    % 提取信道协方差矩阵 RR，然后通过SVD近似宽带波束Space-time block coding
     RR=zeros(N_t,N_t);
     for j=1:num
           RR=RR+1/num*H(:,:,j)'*H(:,:,j);                                 
     end
     [~,~,V]=svd(RR);
     eps = 1e-10;
     F3_temp = V(:,1:N_RF)./ (abs(V(:,1:N_RF))+eps);
     for i_snr = 1:length(SNR)
     for j=1:num
        Fbb3 = F3_temp\F_opt(:,:,j);
        Fbb3 = Fbb3*sqrt(N_s)/norm(F3_temp*Fbb3,'F');
        F3 = F3_temp*Fbb3;
        HH3 = H(:,:,j)*F3;
%         temp3 = temp3+real( log2 ( det(eye(N_s)+HH3'*HH3*SNR )) ) ;
        Rate3(i_snr) = Rate3(i_snr) + real(log2(det(eye(N_s)+HH3'*HH3*SNR(i_snr))));
     end
     end
        
     %% 算法4：本文提出的 DP-AltMin (时相联合优化算法)
     % 联合优化数字预编码、移相器和延时器，逼近 F_opt
    [F_rf, P_DPP, F_BB] = DPALT( F_opt, hc,N_t, N_s, N_RF, num, K, theta, fc, fs);
    % 将优化的延时和移相结合生成最终的模拟预编码矩阵 F_RF
    F_RF = zeros(N_t, N_RF, num);
    for idx_M = 1:num
       F_RF_temp = reshape(F_rf, [P, K*N_RF]); 
       P_temp = P_DPP(:,:,idx_M);
       P_temp = P_temp(:).';
       T_temp = F_RF_temp.*P_temp;
       T_temp = reshape(T_temp, [N_t, N_RF]);
       F_RF(:,:,idx_M) = T_temp;
    end

     for i_snr = 1:length(SNR)
     for j=1:num
        F4_temp = F_RF(:,:,j);
        Fbb4 = F4_temp\F_opt(:,:,j);
        Fbb4 = Fbb4*sqrt(N_s)/norm(F4_temp*Fbb4,'F');
        F4 = F4_temp*Fbb4;
        HH4 = H(:,:,j)*F4;
%         temp2 = temp2+real(log2(det(eye(N_s)+HH2'*HH2*SNR)));
        Rate4(i_snr) = Rate4(i_snr) + real(log2(det(eye(N_s)+HH4'*HH4*SNR(i_snr))));
     end
     end
    %%算法5：TTD-DPP (仅针对径模型的时相混合预编码, 文献[11])
    [F_rf, P_DPP, F_BB] = TTD_DPP( F_opt,N_t, N_s, N_RF, num, K, theta, fc, fs);
    F_RF = zeros(N_t, N_RF, num);
    for idx_M = 1:num
       F_RF_temp = reshape(F_rf, [P, K*N_RF]); 
       P_temp = P_DPP(:,:,idx_M);
       P_temp = P_temp(:).';
       T_temp = F_RF_temp.*P_temp;
       T_temp = reshape(T_temp, [N_t, N_RF]);
       F_RF(:,:,idx_M) = T_temp;
    end

     for i_snr = 1:length(SNR)
     for j=1:num
        F5_temp = F_RF(:,:,j);
        Fbb5 = F5_temp\F_opt(:,:,j);
        Fbb5 = Fbb5*sqrt(N_s)/norm(F5_temp*Fbb5,'F');
        F5 = F5_temp*Fbb5;
        HH5 = H(:,:,j)*F5;
        Rate5(i_snr) = Rate5(i_snr) + real(log2(det(eye(N_s)+HH5'*HH5*SNR(i_snr))));
     end
     end


    %%算法6：MO-Alt (流形优化混合预编码, 文献[13])
     F7_temp = PEAlt(F_opt, N_t, N_RF, N_s);
     for i_snr = 1:length(SNR)
     for j=1:num
        Fbb7 = F7_temp\F_opt(:,:,j);
        Fbb7 = Fbb7*sqrt(N_s)/norm(F7_temp*Fbb7,'F');
        F7 = F7_temp*Fbb7;
        HH7 = H(:,:,j)*F7;
        Rate7(i_snr) = Rate7(i_snr) + real(log2(det(eye(N_s)+HH7'*HH7*SNR(i_snr))));
     end
     end
    toc()
     fprintf('i_loop=%d\n',i);
end

%%4. 取平均并绘图
Rate = Rate/N_iter/num;
Rate1 = Rate1/N_iter/num;
Rate3 = Rate3/N_iter/num;
Rate4 = Rate4/N_iter/num;
Rate5 = Rate5/N_iter/num;
Rate7 = Rate7/N_iter/num;

result = [Rate;Rate4;Rate5;Rate7;Rate3;Rate1;SNR_dB];
%% Figure
N = 6;  
C =linspecer(N);


figure;hold on;grid on; box on;
plot(SNR_dB,Rate,'--', 'color', C(2,:), 'Linewidth',1.8);
plot(SNR_dB,Rate4,'-h','color', C(1,:), 'Linewidth',1.8);
plot(SNR_dB,Rate5,'-v','color', C(3,:), 'Linewidth',1.8);
plot(SNR_dB,Rate7,'-o','color', C(4,:), 'Linewidth',1.8);
plot(SNR_dB,Rate3,'-s','color', C(5,:), 'Linewidth',1.8);
plot(SNR_dB,Rate1,'-d','color', C(6,:), 'Linewidth',1.8);

xlabel('SNR [dB]', 'interpreter', 'latex');
ylabel('Average Rate [bits/s/Hz]', 'interpreter', 'latex');
legend('Fully digital precoding','Proposed DP-AltMin','DPP-TTD [11]', 'MO-AltMin [13]','Space-time block coding [8]', 'Spatially sparse precoding [19]',...
    'interpreter', 'latex');
set(gca,'FontSize',12);
