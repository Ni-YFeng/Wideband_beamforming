function [H,hc,Theta,Beta] = WidebandChannel_Multi(N,M,C,R,fc,fs,tmax,num,sigma_bs, sigma_mt, sigma_t, Theta, Alpha)
%generate the wideband channel with ray-based channel, 1 received antenna
%model, coded by Jingbo tan 2018/11/25
%Transmit antenna number N
%Received antenna number M
%cluster number C
%path number R
%central fre fc
%band fre fs
%max time delay tmax
%subcarrier number num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=zeros(M,N,num+1);

% AODs
Theta = rand(C,1)*pi - pi/2;
sigma=sigma_bs; 
b=sigma/sqrt(2);   
a=rand(C,R)-0.5; 
dTheta = - b*sign(a).*log(1-2*abs(a));

% AOAs
Alpha = rand(C,1)*pi - pi/2;
sigma=sigma_mt; 
b=sigma/sqrt(2);   
a=rand(C,R)-0.5; 
dAlpha = - b*sign(a).*log(1-2*abs(a));


Beta = (randn(C,R)+randn(C,R)*1j)/sqrt(2);

% delays
Delay = rand(C,1)*tmax;
sigma=sigma_t; 
b=sigma/sqrt(2);   
a=rand(C,R)-0.5; 
dDelay = - b*sign(a).*log(1-2*abs(a));

 for k=1:num+1
     if k==num+1
        f=fc;
     else
        f=fc+fs/(num)*(k-1-(num-1)/2);
     end
     for idx_C=1:C
       for idx_R=1:R
       H(:,:,k)=H(:,:,k)+Beta(idx_C,idx_R)*exp(-1j*2*pi*(Delay(idx_C) + dDelay(idx_C, idx_R))*f)*array_response(Alpha(idx_C)+dAlpha(idx_C,idx_R),M,0.5*(3e8)/fc,(3e8)/f)...
           *array_response(Theta(idx_C)+dTheta(idx_C,idx_R),N,0.5*(3e8)/fc,(3e8)/f).';
       end        
     end
 end
H = H.*sqrt(N*M/(C*R));
hc = H(:,:,num+1);
H = H(:,:,1:num);

