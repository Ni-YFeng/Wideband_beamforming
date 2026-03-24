%Test the potential solution for beam squint with multi-stream and multi
%received antennas
%coded by Jingbo Tan 2018/12/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization
L=4;
N=256;
M=4;
fc=3e11;
fs=30e9;
tmax=20e-9;
num=128;
SNR_dB= -20:5:15;
SNR_linear=10.^(SNR_dB/10.);
N_iter=500;
Rate1=zeros(1,length(SNR_dB));
Rate2=zeros(1,length(SNR_dB));
Rate=zeros(1,length(SNR_dB));
Rate3=zeros(1,length(SNR_dB));
Rate4=zeros(1,length(SNR_dB));
Rate5=zeros(1,length(SNR_dB));
H_list=zeros(1,num);
H1_list=zeros(1,num);
H2_list=zeros(1,num);
H3_list=zeros(1,num);
H4_list=zeros(1,num);
%% rate
for i_snr=1:length(SNR_linear)
    SNR=SNR_linear(i_snr);
    temp=0;temp1=0;temp2=0;temp3=0;temp4=0;temp5=0;
    for i =1:N_iter
         [H,hc,Theta, beta] = WidebandChannel_Multi(N,M,L,fc,fs,tmax,num);
         
        %% max
         F_max = zeros(num,N);
         for j=1:num
             [~,~,V]=svd(H(:,:,j));
             F = V(:,1:M);
             F = F/norm(F);
             HH = H(:,:,j)*F;
             temp = temp+real(log2(det(eye(M)+HH*HH'*SNR)));
         end
         
        %% Hybrid near optimal
         R=zeros(N,N);
         for j=1:num
               R=R+1/num*H(:,:,j)'*H(:,:,j);
         end
         [~,~,V]=svd(R);
         F1_temp = V(:,1:M);
         for j=1:num
            H1_temp = H(:,:,j)*F1_temp;
            [~,~,Vbb] = svd(H1_temp);
            Fbb = Vbb(:,1:M);
            F1 = F1_temp*Fbb/norm(F1_temp*Fbb);
            HH1 = H(:,:,j)*F1;
            temp1 = temp1+real(log2(det(eye(M)+HH1*HH1'*SNR)));
         end 
         
        %% Wide beam
         for j=1:M
              f1=fc+fs/(num)*(-(num-1)/2);
              f2=fc+fs/(num)*((num-1)/2);
              width = abs(f2/fc*sin(Theta(j))-f1/fc*sin(Theta(j)));
              F2_temp(:,j) = GenerateWideBeam(N,sin(Theta(j)),width);
         end
         for j=1:num
             H2_temp = H(:,:,j)*F2_temp;
             [~,~,Vbb] = svd(H2_temp);
             Fbb = Vbb(:,1:M);
             F2 = F2_temp*Fbb/norm(F2_temp*Fbb);
             HH2 = H(:,:,j)*F2;
             temp2 = temp2+real(log2(det(eye(M)+HH2*HH2'*SNR)));
         end 
         
        %% Partial TTD
         Sub_num = 16;
         for j=1:M
             Phase(j,:) =[0:N-1]*pi*sin(Theta(j));
             Delay(j) = N/Sub_num*sin(Theta(j))/2;
         end
         for j=1:num
            f=fc+fs/(num)*(j-1-(num-1)/2);
            for k=1:M
                beta = 2*(f/fc-1)*Delay(k);
                F3_temp(k,:) = GenerateComBeam(N,Sub_num,beta,sin(Theta(k)));
            end
            H3_temp = H(:,:,j)*F3_temp.';
            [~,~,Vbb] = svd(H3_temp);
            Fbb = Vbb(:,1:M);
            F3 = F3_temp.'*Fbb/norm(F3_temp.'*Fbb);
            HH3 = H(:,:,j)*F3;
            temp3 = temp3+real(log2(det(eye(M)+HH3*HH3'*SNR)));
         end 
        %% Traditional
         theta_list = -1:0.01:1;
         for j = 1:length(theta_list)
            At(j,:)=1/sqrt(N)*exp(-1j*[0:N-1]*pi*theta_list(j));
         end
         [F4_temp,F_BB]=spatially_sparse_precoding(M,M,hc,At.');
         for j=1:num
            H4_temp = H(:,:,j)*F4_temp;
            [~,~,Vbb] = svd(H4_temp);
            Fbb = Vbb(:,1:M);
            F4 = F4_temp*Fbb/norm(F4_temp*Fbb);
            HH4 = H(:,:,j)*F4;
            temp4 = temp4+real(log2(det(eye(M)+HH4*HH4'*SNR)));
         end 
         
         fprintf('i_SNR=%d\n,i_loop=%d\n',i_snr,i);
    end
    Rate(i_snr)=real(temp/(N_iter*num));
    Rate1(i_snr)=real(temp1/(N_iter*num));
    Rate2(i_snr)=real(temp2/(N_iter*num));
    Rate3(i_snr)=real(temp3/(N_iter*num));
    Rate4(i_snr)=real(temp4/(N_iter*num));
    Rate5(i_snr)=real(temp5/(N_iter*num));
end
%% Figure

figure;
plot(SNR_dB,Rate,'--','Linewidth',1.2,'Markersize',5);
hold on;
%plot(SNR_dB,Rate5,'-<','Linewidth',1.2);
%hold on;
plot(SNR_dB,Rate3,'-o','Linewidth',1.2,'Markersize',5);
%hold on;
%plot(SNR_dB,Rate5,'-<','Linewidth',1.5);
hold on;
plot(SNR_dB,Rate1,'->','Linewidth',1.2,'Markersize',5);
hold on;
plot(SNR_dB,Rate2,'-v','Linewidth',1.2,'Markersize',5);
hold on;
plot(SNR_dB,Rate4,'-d','Linewidth',1.2,'Markersize',5);


grid on;
l1=xlabel('SNR (dB)');
l2=ylabel('Acheivabel rate per subcarrier (bit/s/Hz)');
l3=legend('Optimal precoding', 'Proposed TTD-DPP','Optimizaiton based hybrid precoding [16]','Wide beam hybrid precoding [19]','Spartially sparse precoding [14]');
set(l1,'FontSize',12);
set(l2,'FontSize',12);
set(l3,'FontSize',10);