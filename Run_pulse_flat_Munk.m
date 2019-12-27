close all;
clear variables;
clc;

set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial'); 


% source function

dt = 10^(-4);
Tmax = 0.5;
tt = 0:dt:Tmax;

alph = 400;
bet = 500;
aa = (bet^2+alph^2)/alph;
phi = pi/2;

ft_s = aa*exp(-bet*tt).*(-bet*sin(alph*tt + phi) + alph*cos(alph*tt + phi) + bet*sin(phi) )  ;
% 
% figure;
% plot(tt,ft_s,'linewidth',1.5);
% grid on;
% legend('dg(t)/dt');
% xlabel('t, s');
% xlim([0 0.1]);

% 
% gt_s = aa*exp(-bet*tt).*(sin(alph*tt+phi) - sin(phi)) ;
% 
% figure;
% plot(tt,gt_s,'linewidth',1.5,'color','red');
% grid on;
% legend('g(t)');
% xlabel('t, s');


% source-receiver geometry

zs = 2000;
zr = 2000;
R = [5000 25000 50000 75000 100000];

% medium

MP.LayersData = [[0 1500 1500 1 1 0 0]; [4500 1500 1600 1 1 0 0] ];
z = 0:2:5000;
zbar = 2*(z-1300)/1300;
eps = 0.00737;
c = 1500*(1 + eps*(zbar - 1 + exp(-zbar)));

figure;
plot(c,z/1000,'linewidth',1.5);
grid on;
xlabel('c, m/s');
ylabel('z, km');


MP.HydrologyData = [z.' c.'];


[t_r pt_r] = PulseSignalAssembleNoR(  ft_s, dt, zs, zr, MP, R  );


for ii = 1:length(R)
    
    dlmwrite(['pulses/' int2str(R(ii)) '_ac_modes_bb.txt'],[t_r(:,ii) pt_r(:,ii)],'delimiter','\t','precision',8);
    
end;



