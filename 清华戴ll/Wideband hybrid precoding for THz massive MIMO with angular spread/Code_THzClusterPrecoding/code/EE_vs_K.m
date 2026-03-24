clear all;
%% initialization
C= 2;
R = 10;
N_t = 256;  %antenna number
N_s = 2;
N_r = 2;  
N_RF = 2;

fc=10e10;
fs=0.1*fc;
tmax=25e-9;
sigma_t = 0.06e-9;
num = 32;   %frequence number
SNR_dB= 10;
SNR = 10.^(SNR_dB/10.);

N_iter = 100;
sigma_bs = 3 * pi / 180;
sigma_mt = 3 * pi / 180;

% 遍历不同的延时器数量 K
K_list = [1,2,4,8,16,32];
% P = N_t/K;

 theta_list = linspace(-1,1,512)*pi/2;
 At = zeros(length(theta_list), N_t);
 for j = 1:length(theta_list)
    At(j,:)=1/sqrt(N_t)*exp(-1j*[0:N_t-1]*pi*sin(theta_list(j)));
 end
 
Rate=zeros(1,length(K_list));
Rate1=zeros(1,length(K_list));
Rate2=zeros(1,length(K_list));
Rate3=zeros(1,length(K_list));
Rate4=zeros(1,length(K_list));
Rate5=zeros(1,length(K_list));
Rate6=zeros(1,length(K_list));
Rate7=zeros(1,length(K_list));


temp0=0;temp1=0;temp2=0;temp3=0;temp4=0;temp5=0; temp6 = 0;temp7 = 0;
parfor i1 = 1:length(K_list)
            tic()
    K = K_list(i1);
    P = N_t/K;% 更新每个延时器连接的移相器数量
    %% 定义不同架构的系统功耗模型 (单位: W)
    % 1. 全数字预编码功耗 (FD): 发射功耗 + 基带功耗 + N_t个射频链功耗
    E_FD = (1000 + 300 + N_t*300)/1000;
    % 2. 提出的时相混合预编码功耗 (DP): 
    % 发射 + 基带 + N_RF个射频链 + N_RF*K个延时器(100mW) + N_RF*N_t个移相器(40mW
    E_DP = (1000 + 300 + N_RF*300 + N_RF*K*100 + N_RF*N_t*40)/1000;
    % 3. 传统混合预编码功耗 (HP): 
    % 发射 + 基带 + N_RF个射频链 + N_RF*N_t个移相器(无延时器)
    E_HP = (1000 + 300 + N_RF*300 + N_RF*N_t*40)/1000;
    for i2 =1:N_iter

           % 生成信道和F_opt
         [H,hc,Theta, beta] = WidebandChannel_Multi(N_t,N_r,C,R,fc,fs,tmax,num,sigma_bs, sigma_mt,sigma_t);
         F_opt = zeros(N_t, N_s, num);
         for j=1:num
             [~,~,V]=svd(H(:,:,j));
             F_opt(:,:,j) = V(:,1:N_s);
         end
         [~,~,V]=svd(hc);
        V = V(:,1:N_s);
        F_optc = V;
        %% 计算速率并直接除以对应的功耗，得到能量效率
        %% optimal

         F_max = zeros(num,N_t);
         
         for j=1:num
             F = F_opt(:,:,j);
             HH = H(:,:,j)*F;
             Rate(i1) = Rate(i1) + real(log2(det(eye(N_s)+HH'*HH*SNR)))/E_FD;
         end


        %% OMP
         [F1_temp,F_BB,th]=spatially_sparse_precoding(N_RF,N_s,F_optc,At.');
         th = -theta_list(th);
         Theta = th;
         for j=1:num
            Fbb1 = F1_temp\F_opt(:,:,j);
            Fbb1 = Fbb1*sqrt(N_s)/norm(F1_temp*Fbb1,'F');
            F1 = F1_temp*Fbb1;
            HH1 = H(:,:,j)*F1;
            Rate1(i1) = Rate1(i1)+real(log2(det(eye(N_s)+HH1'*HH1*SNR)))/E_HP;
         end
         
        %% wideband
         RR=zeros(N_t,N_t);
         for j=1:num
               RR=RR+1/num*H(:,:,j)'*H(:,:,j);                                 
         end
         [~,~,V]=svd(RR);
         eps = 1e-10;
         F3_temp = V(:,1:N_RF)./ (abs(V(:,1:N_RF))+eps);
         for j=1:num
            Fbb3 = F3_temp\F_opt(:,:,j);
            Fbb3 = Fbb3*sqrt(N_s)/norm(F3_temp*Fbb3,'F');
            F3 = F3_temp*Fbb3;
            HH3 = H(:,:,j)*F3;
    %         temp3 = temp3+real( log2 ( det(eye(N_s)+HH3'*HH3*SNR )) ) ;
            Rate3(i1) = Rate3(i1) + real(log2(det(eye(N_s)+HH3'*HH3*SNR)))/E_HP;
         end
         

        %% TTD-ALT
        [F_rf, P_DPP, F_BB] = DPALT( F_opt, hc,N_t, N_s, N_RF, num, K, Theta, fc, fs);
        F_RF = zeros(N_t, N_RF, num);
        for idx_M = 1:num
           F_RF_temp = reshape(F_rf, [P, K*N_RF]); 
           P_temp = P_DPP(:,:,idx_M);
           P_temp = P_temp(:).';
           T_temp = F_RF_temp.*P_temp;
           T_temp = reshape(T_temp, [N_t, N_RF]);
           F_RF(:,:,idx_M) = T_temp;
        end
 
         for j=1:num
            F4_temp = F_RF(:,:,j);
            Fbb4 = F4_temp\F_opt(:,:,j);
            Fbb4 = Fbb4*sqrt(N_s)/norm(F4_temp*Fbb4,'F');
            F4 = F4_temp*Fbb4;
            HH4 = H(:,:,j)*F4;
            Rate4(i1) = Rate4(i1)+real( log2 ( det(eye(N_s)+HH4'*HH4*SNR)))/E_DP;
         end
         
       %% TTD
        [F_rf, P_DPP, F_BB] = TTD_DPP( F_opt,N_t, N_s, N_RF, num, K, Theta, fc, fs);
        F_RF = zeros(N_t, N_RF, num);
        for idx_M = 1:num
           F_RF_temp = reshape(F_rf, [P, K*N_RF]); 
           P_temp = P_DPP(:,:,idx_M);
           P_temp = P_temp(:).';
           T_temp = F_RF_temp.*P_temp;
           T_temp = reshape(T_temp, [N_t, N_RF]);
           F_RF(:,:,idx_M) = T_temp;
        end
 
         for j=1:num
            F5_temp = F_RF(:,:,j);
            Fbb5 = F5_temp\F_opt(:,:,j);
            Fbb5 = Fbb5*sqrt(N_s)/norm(F5_temp*Fbb5,'F');
            F5 = F5_temp*Fbb5;
            HH5 = H(:,:,j)*F5;
            Rate5(i1) = Rate5(i1)+real( log2 ( det(eye(N_s)+HH5'*HH5*SNR)))/E_DP;
         end
         
        %% MOALT
         F6_temp = PEAlt(F_opt, N_t, N_RF, N_s);
         for j=1:num
            Fbb6 = F6_temp\F_opt(:,:,j);
            Fbb6 = Fbb6*sqrt(N_s)/norm(F6_temp*Fbb6,'F');
            F6 = F6_temp*Fbb6;
            HH6 = H(:,:,j)*F6;
            Rate6(i1) = Rate6(i1)+real( log2 ( det(eye(N_s)+HH6'*HH6*SNR )) )/E_HP ;
         end
        toc()
         fprintf('[%d/%d] i_loop=%d\n',i1, length(K_list), i2);
    end
end

Rate = Rate/N_iter/num;
Rate1 = Rate1/N_iter/num;
Rate3 = Rate3/N_iter/num;
Rate4 = Rate4/N_iter/num;
Rate5 = Rate5/N_iter/num;
Rate6 = Rate6/N_iter/num;

% result = [Rate;Rate4;Rate5;K_list];
%% Figure
N = 6;  
C =linspecer(N);


figure;hold on;grid on; box on;
plot(K_list, Rate,'--', 'color', C(2,:), 'Linewidth',1.8);
plot(K_list, Rate4,'-h','color', C(1,:), 'Linewidth',1.8, 'markerfacecolor', 'w');
plot(K_list, Rate5,'-v','color', C(3,:), 'Linewidth',1.8, 'markerfacecolor', 'w');
plot(K_list, Rate6,'-o','color', C(4,:), 'Linewidth',1.8, 'markerfacecolor', 'w');
plot(K_list, Rate3,'-s','color', C(5,:), 'Linewidth',1.8, 'markerfacecolor', 'w');
plot(K_list, Rate1,'-d','color', C(6,:), 'Linewidth',1.8, 'markerfacecolor', 'w');

xlabel('$K$',    'interpreter', 'latex');
ylabel('Energy efficiency [bps/Hz/W]',    'interpreter', 'latex');
legend('Fully digital precoding','Proposed DP-AltMin','DPP-TTD [11]', 'MO-AltMin [13]','Space-time block coding [8]', 'Spatially sparse precoding [19]',...
    'interpreter', 'latex');
set(gca,'FontSize',12);
xlim([1, 32])