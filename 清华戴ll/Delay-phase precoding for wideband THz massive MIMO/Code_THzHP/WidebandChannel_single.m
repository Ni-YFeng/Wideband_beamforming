function [H,hc,Theta,Beta] = WidebandChannel_single(N,L,fc,fs,tmax,num)
%generate the wideband channel with ray-based channel, 1 received antenna
%model, coded by Jingbo tan 2018/11/25
%antenna number N
%path number L
%central fre fc
%band fre fs
%max time delay tmax
%subcarrier number num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=zeros(N,num+1);
%Theta = rand(L,1)*pi-pi/2;
Theta = asin(0.5);
%Beta = randn(L,1)+randn(L,1)*1j;
Beta = ones(L);
Delay = rand(L,1)*tmax;
 for k=1:num+1
     if k==num+1
        f=fc;
     else
        f=fc+fs/(num)*(k-1-(num-1)/2);
     end
     for j=1:L        
       H(:,k)=H(:,k)+Beta(j)*exp(-1j*2*pi*Delay(j)*f)*array_respones(Theta(j),N,0.5*(3e8)/fc,(3e8)/f);
     end        
 end
H = H.*sqrt(N/L);
hc = H(:,num+1);
H = H(:,1:num);

