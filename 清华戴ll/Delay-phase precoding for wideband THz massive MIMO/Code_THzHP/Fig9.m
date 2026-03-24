%Draw the array gain comparison of narrow and wide band
H_n = WidebandChannel_single(256,1,3e11,0.3e9,0,128);
H_n1 = WidebandChannel_single(256,1,3e11,3e9,0,128);
H_w = WidebandChannel_single(256,1,3e11,30e9,0,128);
H_w1 = WidebandChannel_single(64,1,28e9,2.8e9,0,128);
fc = 3e11;
f=1/sqrt(256)*exp(1j*[0:256-1]*pi*0.5);
f1=1/sqrt(64)*exp(1j*[0:64-1]*pi*0.5);
Array_gain_sub=zeros(1,128);for i=1:128
    fm=3e11+30e9/128*(i-1-(128-1)/2);
    beta = 2*(fm/fc-1)*256/16*sin(pi/6)/2;
    f2=GenerateComBeam(256,16,beta,sin(pi/6));
    Array_gain_sub(i)=abs(f2*H_w(:,i));
end
Array_gain_n = abs(f*H_n)./max(abs(f*H_n));
Array_gain_n1 = abs(f*H_n1)./max(abs(f*H_n1));
Array_gain_w = abs(f*H_w)./max(abs(f*H_w));
Array_gain_w1 = abs(f1*H_w1)./max(abs(f1*H_w1));
Array_gain_sub = Array_gain_sub./max(Array_gain_sub);
figure;
%subplot(2,1,2)
plot([1:1:128],Array_gain_n, '--','Linewidth', 1);
hold on;
plot([1:1:128],Array_gain_n1, '-.','Linewidth', 1);
hold on;
plot([1:1:128],Array_gain_w,'.-','Linewidth', 1,'MarkerSize',10);
%hold on;
%plot([1:1:128],Array_gain_w1 ,'Linewidth', 1.2);

hold on;
plot([1:1:128],Array_gain_sub, 'Linewidth', 1);
grid on;
ylim([0,1.1]);
xlim([1,128]);
l2=xlabel('Subcarrier (m)');
l3=ylabel('Normalized array gain');
l1=legend('Hybrid Precoding, B=0.1 GHz','Hybrid Precoding, B=1 GHz','Hybrid Precoding, B=10 GHz','Proposed DPP, B=10 GHz');
set(l1,'FontSize',10);
set(l2,'FontSize',12);
set(l3,'FontSize',12);
