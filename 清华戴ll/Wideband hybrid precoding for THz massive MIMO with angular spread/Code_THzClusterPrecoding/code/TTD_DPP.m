 function [F_RF, P_DPP, F_BB] = TTD_DPP(F_opt, N_t, N_s, N_RF, M, K, Theta, fc, fs)

 
P = N_t/K;


 for k=1:N_RF
     F_RF(:,k) = exp(1i*pi*(0:N_t-1)*sin(Theta(k))).';
 end
phi = angle(F_RF);

%% ¼ÆËãF_BB
F_BB = randn(N_RF,N_s,M);
 for idx_M=1:M
    F_BB(:,:,idx_M) = F_RF\F_opt(:,:,idx_M);
 end

eta = zeros(M);

t = rand(K,N_RF);
for idx = 1:N_RF
    theta_c = sin(Theta(idx));
    s = -P*theta_c/2;
    if theta_c > 0
        t(:,idx) = (K-1)*abs(s) + (1:K)'*ceil(s); 
    else
        t(:,idx) = (1:K)'*ceil(s); 
    end
end

P_DPP = zeros(K,N_RF,M);
for m = 1:M
   fm = fc+fs/(M)*(m-1-(M-1)/2);
   eta(m) = fm/fc;
   P_DPP(:,:,m) = exp(-1j*2*pi*eta(m)*t); 
end


 


 end
 


