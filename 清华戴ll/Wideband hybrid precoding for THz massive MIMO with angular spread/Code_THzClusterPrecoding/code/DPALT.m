 function [F_RF, P_DPP, F_BB] = DPALT(F_opt, hc, N_t, N_s, N_RF, M, K, Theta, fc, fs)

 
P = N_t/K;
[~,~,V]=svd(hc);
V = V(:,1:N_s);

%% 初始化F_RF
 for k=1:N_RF
     F_RF(:,k) = exp(1i*pi*(0:N_t-1)*sin(Theta(k))).';
 end
phi = angle(F_RF);

%% 初始化F_BB
F_BB = randn(N_RF,N_s,M);
 for idx_M=1:M
    F_BB(:,:,idx_M) = F_RF\F_opt(:,:,idx_M);
 end

%% 初始化P_DPP


eta = zeros(M,1);
for m =1:M
   fm = fc+fs/(M)*(m-1-(M-1)/2);
   eta(m) = fm/fc;
end
label = zeros(K, M, N_RF);
R = (1:K)';
for idx = 1:N_RF
    theta_c = sin(Theta(idx));
    b = sin(P*pi/2*(1-eta')*theta_c)./sin(pi/2*(1-eta')*theta_c);
    a = pi*theta_c*((R-0.5)*P-0.5)*(1-eta');
    label(:,:,idx) = exp(1j*a).*b;
end
t2 = zeros(K, N_RF);
S = 1024;
tcode = linspace(-N_t/2, N_t/2, S);
Tcode = zeros(M, length(tcode));
for idx = 1:length(tcode)
   Tcode(:,idx) = exp(-1j*2*pi*tcode(idx)*eta); 
end

Value = [];
Index = [];
for idx1 = 1:N_RF
   for idx2 = 1:K
      Q = label(idx2,:,idx1);
      pha = real(Q*Tcode);
      [value,index]=max(pha);
      Value = [Value, value];
      Index = [Index, index];
      t2(idx2, idx1) = tcode(index);
   end
end
P_DPP = zeros(K,N_RF,M);
for idx = 1:N_RF
   if min(t2(:, idx)) <0 
        t2(:,idx) = t2(:,idx) - min(t2(:,idx));
   end 
end
t = t2;

for idx = 1:M
   P_DPP(:,:,idx) = exp(-1j*2*pi*eta(idx)*t); 
end

%% LALT
N_iter = 2;
NN = 10;
for nn = 1:NN
    
    %% Update Delay
    t3 = zeros(K,N_RF);
    label = zeros(K,N_RF, M);
    for idx_M = 1:M
       temp = (F_opt(:,:,idx_M)*pinv(F_BB(:,:,idx_M))).*conj(F_RF);
       temp = reshape(temp, [P, K*N_RF]);
       temp = sum(temp,1);
       temp = reshape(temp, [K,N_RF]);
       label(:,:,idx_M) = temp;
    end
    for idx1 = 1:K
       for idx2 = 1:N_RF
          Q = label(idx1,idx2,:);
          Q = squeeze(Q);
          pha = real(Q'*Tcode);
          [value,index]=max(pha);
          t3(idx1, idx2) = tcode(index);
       end
    end
    P_DPP = zeros(K,N_RF,M);
    for idx = 1:N_RF
       if min(t3(:, idx)) <0 
            t3(:,idx) = t3(:,idx) - min(t3(:,idx));
       end 
    end
    for idx = 1:M
       P_DPP(:,:,idx) = exp(-1j*2*pi*eta(idx)*t3); 
    end
    

for n_iter = 1:N_iter
    
        %% update digital precoder
        for idx_M = 1:M
           P_temp = kron(P_DPP(:,:,idx_M), ones(P,1));
           T_temp = F_RF.*P_temp;
            F_BB(:,:,idx_M) =T_temp\F_opt(:,:,idx_M);
        end
        
        %% update analog phase shift

           F_temp = zeros(size(F_RF));
           for idx_M = 1:M
               P_temp = kron(P_DPP(:,:,idx_M), ones(P,1));
               F_temp =  F_temp + norm(F_BB(:,:,idx_M),'F')^2*(F_opt(:,:,idx_M)*pinv(F_BB(:,:,idx_M))).*conj(P_temp);
           end
           F_RF = F_temp./abs(F_temp);
            

end
 


end
 end





