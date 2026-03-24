clear all;
%% initialization
C = 1;
R=  10;
N_t=256;  %antenna number
N_s= 1;  
N_r = 1;
N_RF = 1;

fc=1e11;
fs=0.1*fc;
tmax=25e-9;
sigma_t = 0.06e-9;
num=64;   %frequence number
SNR_dB= 10;
SNR=10.^(SNR_dB/10.);
N_iter= 300;
sigma_bs = 3*pi/180;
sigma_mt = 3*pi/180;
K = 16;
P = N_t/K;
theta_list = linspace(-1,1,2048)*pi/2;
At = zeros(length(theta_list), N_t);
for j = 1:length(theta_list)
At(j,:)=1/sqrt(N_t)*exp(-1j*[0:N_t-1]*pi*sin(theta_list(j)));
end
X = linspace(0*pi/180, 5*pi/180, 8);
X2 = linspace(0, 0.06e-9, 8);
% X = 0
% X = -20:5:15;
Rate=zeros(1,length(X));
Rate1=zeros(1,length(X));
Rate2=zeros(1,length(X));
Rate3=zeros(1,length(X));
Rate4=zeros(1,length(X));
Rate5=zeros(1,length(X));
Rate6=zeros(1,length(X));
%% rate
parfor idx_X=1:length(X)
    fprintf('idx=%d\n',idx_X);
    sigma_bs = X(idx_X);
    sigma_mt = sigma_bs;
%     SNR = 10.^(X(idx_X)/10);
    sigma_t = X2(idx_X);
    temp0=0;temp1=0;temp2=0;temp3=0;temp4=0;temp6 = 0;temp5=0;
    for i =1:N_iter
         tic()
%          disp('channel')
%          rng(i);
         [H,hc,Theta, beta] = WidebandChannel_Multi(N_t,N_r,C,R,fc,fs,tmax,num,sigma_bs, sigma_mt,sigma_t);
         F_opt = zeros(N_t, N_s, num);
         for j=1:num
             [~,~,V]=svd(H(:,:,j));
             F_opt(:,:,j) = V(:,1:N_s);
         end
         [~,~,V]=svd(hc);
        V = V(:,1:N_s);
        F_optc = V;
        %% optimal

         F_max = zeros(num,N_t);
         
         for j=1:num
             F = F_opt(:,:,j);
             HH = H(:,:,j)*F;
             temp0 = temp0+real(log2(det(eye(N_s)+HH'*HH*SNR)));
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
            temp1 = temp1+real(log2(det(eye(N_s)+HH1'*HH1*SNR)));
         end
         

         

        %% DPALT
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
            temp4 = temp4+real( log2 ( det(eye(N_s)+HH4'*HH4*SNR)));
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
            temp5 = temp5+real( log2 ( det(eye(N_s)+HH5'*HH5*SNR)));
         end
         
        %% MOALT
         F6_temp = PEAlt(F_opt, N_t, N_RF, N_s);
         for j=1:num
            Fbb6 = F6_temp\F_opt(:,:,j);
            Fbb6 = Fbb6*sqrt(N_s)/norm(F6_temp*Fbb6,'F');
            F6 = F6_temp*Fbb6;
            HH6 = H(:,:,j)*F6;
            temp6 = temp6+real( log2 ( det(eye(N_s)+HH6'*HH6*SNR )) ) ;
         end
         
         
        toc()
         fprintf('idx=%d, i_loop=%d\n',idx_X, i);
    end
    Rate(idx_X)=real(temp0/(N_iter*num));
    Rate1(idx_X)=real(temp1/(N_iter*num));
    Rate2(idx_X)=real(temp2/(N_iter*num));
    Rate4(idx_X)=real(temp4/(N_iter*num));
    Rate5(idx_X)=real(temp5/(N_iter*num));
    Rate6(idx_X)=real(temp6/(N_iter*num));
end
%% Figure
N = 6;  
C =linspecer(N);
X1 = X*180/pi;
result = [Rate; Rate4;Rate5;Rate6;Rate1;X1];
figure;hold on;grid on;box on;
plot(X1,Rate,'--','color', C(2,:),'Linewidth',1.6);
plot(X1,Rate4,'-v','color', C(1,:),'Linewidth',1.6);
plot(X1,Rate5,'-<','color', C(3,:),'Linewidth',1.6);
plot(X1,Rate6,'-s','color', C(4,:),'Linewidth',1.6);
plot(X1,Rate1,'-o','color', C(6,:),'Linewidth',1.6);

xlabel('Angular Spread [degree]', 'interpreter', 'latex');
ylabel('Average Rate [bits/s/Hz]', 'interpreter', 'latex');
legend('Fully digital precoding', 'Proposed DP-AltMin', 'DPP-TTD [11]', 'MO-AltMin [13]', 'Spatially sparse precoding [19]', 'interpreter', 'latex');
set(gca,'FontSize',12);