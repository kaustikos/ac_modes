close all;
clear variables;
clc;

xx = dlmread('pulses/ft.txt');
dt = xx(2,1)-xx(1,1);
ft_s = xx(:,2);
zs = 50;

zr = 10;
R = [5000 6500 8000];

MP = [[0 1500 1500 1 1]; [90 1500 1800 1 1.8]];



[t_r pt_r] = PulseSignalAssemble(  ft_s, dt, zs, zr, MP, R  );


for ii = 1:length(R)
    
    dlmwrite(['pulses/' int2str(R(ii)) '_ac_modes.txt'],[t_r(:,ii) pt_r(:,ii)],'delimiter','\t','precision',8);
    
end;