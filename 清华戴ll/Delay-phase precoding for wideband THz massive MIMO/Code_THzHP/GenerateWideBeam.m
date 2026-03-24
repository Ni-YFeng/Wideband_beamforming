function [b] = GenerateWideBeam(N,phi,w)
%Generate wide with specific center and width
%coded by Jingbo 2018/11/26
%antenna number N
%center phi
%beam width w
width=w*N/(N-1);
p=exp(1j*pi*[1:N]*phi);
q=exp(-1j*pi*[1:N].*([1:N]-N+1)*width/N);
b=p.*q;
end

