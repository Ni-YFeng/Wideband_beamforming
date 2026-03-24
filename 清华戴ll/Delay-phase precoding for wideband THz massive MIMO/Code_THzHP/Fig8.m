%Draw the Physical direction of Time Phase beamforming
theta_list = -1:0.0002:1;
f_list = [1;128/2;128];
Angle_table1 = zeros(length(theta_list),256);
Angle_table2 = zeros(length(theta_list),256);
Angle_table3 = zeros(length(theta_list),256);
for i = 1:length(theta_list)
    Angle_table1(i,:)=1/sqrt(256)*exp(1j*[0:256-1]*pi*285/300*theta_list(i));
    Angle_table2(i,:)=1/sqrt(256)*exp(1j*[0:256-1]*pi*300/300*theta_list(i));
    Angle_table3(i,:)=1/sqrt(256)*exp(1j*[0:256-1]*pi*315/300*theta_list(i));
end
Delay = 256/16*0.5/2;
beta1 = 2*(285/300-1)*Delay;
beta2 = 2*(300/300-1)*Delay;
beta3 = 2*(315/300-1)*Delay;
a1 = GenerateComBeam(256,16,beta1,0.5);
a2 = GenerateComBeam(256,16,beta2,0.5);
a3 = GenerateComBeam(256,16,beta3,0.5);
Array_gain_1 = abs(Angle_table1*a1')./max(abs(Angle_table2*a2'));
Array_gain_c = abs(Angle_table2*a2')./max(abs(Angle_table2*a2'));
Array_gain_M = abs(Angle_table3*a3')./max(abs(Angle_table2*a2'));

figure;
subplot(1,2,1);
plot(theta_list,Array_gain_1,'--');
hold on;
plot(theta_list,Array_gain_c);
hold on;
plot(theta_list,Array_gain_M,'-.');
grid on;
xlim([0.4,0.6]);
l1=xlabel({'$\theta$'},'Interpreter','latex');
l2=ylabel('Normalized array gain');
legend({'$f_{1}$','$f_{c}$','$f_{M}$'},'Interpreter','latex');
set(l1,'FontSize',12);
set(l2,'FontSize',12);
