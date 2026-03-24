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

N_iter = 500;
sigma_bs = 5 * pi / 180;
sigma_mt = 5 * pi / 180;

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
    K = K_list(i1);
    P = N_t/K;
    for i2 =1:N_iter
        tic()
         [H,hc, theta] = WidebandChannel_Multi(N_t,N_r,C,R,fc,fs,tmax,num,sigma_bs, sigma_mt,sigma_t);
         F_opt = zeros(N_t, N_s, num);
         for j=1:num
             [~,~,V]=svd(H(:,:,j));
             F_opt(:,:,j) = V(:,1:N_s);
         end
         [~,~,V]=svd(hc);
        V = V(:,1:N_s);
        F_optc = V;
        
        [~,~,th]=spatially_sparse_precoding(N_RF,N_s,F_optc,At.');
        theta = -theta_list(th);
        %% optimal
         F_max = zeros(num,N_t);


         for j=1:num
             F = F_opt(:,:,j);
             HH = H(:,:,j)*F;
             Rate(i1) = Rate(i1) + real(log2(det(eye(N_s)+HH'*HH*SNR)));
         end



       %% DPALT
        [F_rf, P_DPP, F_BB] = DPALT( F_opt, hc,N_t, N_s, N_RF, num, K, theta, fc, fs);
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
    %         temp2 = temp2+real(log2(det(eye(N_s)+HH2'*HH2*SNR)));
            Rate4(i1) = Rate4(i1) + real(log2(det(eye(N_s)+HH4'*HH4*SNR)));
         end

        %% TTD
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


         for j=1:num
            F5_temp = F_RF(:,:,j);
            Fbb5 = F5_temp\F_opt(:,:,j);
            Fbb5 = Fbb5*sqrt(N_s)/norm(F5_temp*Fbb5,'F');
            F5 = F5_temp*Fbb5;
            HH5 = H(:,:,j)*F5;
            Rate5(i1) = Rate5(i1) + real(log2(det(eye(N_s)+HH5'*HH5*SNR)));
         end

        toc()
         fprintf('[%d/%d] i_loop=%d\n',i1, length(K_list), i2);
    end
end

Rate = Rate/N_iter/num;
Rate4 = Rate4/N_iter/num;
Rate5 = Rate5/N_iter/num;

result = [Rate;Rate4;Rate5;K_list];
%% Figure
N = 6;  
C =linspecer(N);


figure;hold on;grid on; box on;
plot(K_list,Rate,'--', 'color', C(2,:), 'Linewidth',1.8);
plot(K_list,Rate4,'-o','color', C(1,:), 'Linewidth',1.8, 'markerfacecolor', 'w');
plot(K_list,Rate5,'-v','color', C(3,:), 'Linewidth',1.8, 'markerfacecolor', 'w');

xlabel('$K$', 'interpreter', 'latex');
ylabel('Average rate [bps/Hz]', 'interpreter', 'latex');
legend('Fully digital precoding','Proposed DP-AltMin','DPP-TTD [11]', 'interpreter', 'latex');
set(gca,'FontSize',12);
xlim([1, 32])
ylim([8, 22])