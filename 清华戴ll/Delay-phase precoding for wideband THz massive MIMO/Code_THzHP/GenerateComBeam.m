function [w] = GenerateComBeam(N,K,beta, alpha)
%Generate the combined beam
%coded by Jingbo 2018/11/27
%N antennas number
%K sub-beam size
%beta sub-beam angle
%alpha specific angle
B1 = exp(1j*pi*[0:N-1]*alpha);
B2 = exp(1j*pi*[0:K-1]*beta);
M = N/K;
w = zeros(1,N);
for i=1:K
    w((i-1)*M+1:i*M)=B1((i-1)*M+1:i*M)*B2(i);
end

