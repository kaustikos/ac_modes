close all;
clear variables;
clc;


set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial');

cw = 1500;
cb = 2000;

rhow = 1;
rhob = 2;

betab = 0.5;

h0 = 50;

freq = 100;

omeg = 2*pi*freq;

dz0 = 0.125;
opts.Tgr = 3;
opts.Ngr = 3;
opts.nmod = 20;
opts.Hb = 200;
opts.BotBC = 'D';

dhs = 0;

kh = [];

% c(z) = c_0 – (\Delta c) tanh( (z-z_0)/\sigma),

%c0 = 1487.8 m/s, z_0 = 32 m, \sigma = 13 m, \Delta c = 33.17 m/s.

c0 = 1490;
z0 = 25;
dc = 30;
sigma = 10;

zssp = 0:2:100;
cssp = c0 - dc*tanh( (zssp-z0)/sigma );

MP.HydrologyData = [zssp.' cssp.'];

figure;
plot(cssp,zssp);
set(gca,'Ydir','reverse');


MP.LayersData = [[0    cw  cw  rhow    rhow   0    0];
                          [h0   cw  cb  rhow    rhob   0    betab]];


[krs, wmode] = ac_modesr(dz0,MP,freq,opts);

z = dz0*(0:size(wmode,1)-1);

figure;
hold all;
for ii = 1:opts.nmod
    plot(z,wmode(:,ii));
end;

kj_im = ModesAttCoeffs(dz0,freq,krs,wmode,MP);
mgv = ModesGroupVelocities(z,freq,krs,wmode,MP);

figure;
plot(krs+1i*kj_im,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);


dlmwrite(['case_2/kj_f_' int2str(freq) 'Hz_h_50m.txt'], (krs+kj_im).','delimiter','\t','precision',10);
dlmwrite(['case_2/vgj_f_' int2str(freq) 'Hz_h_50m.txt'], (mgv).','delimiter','\t','precision',10);
dlmwrite(['case_2/phij_f_' int2str(freq) 'Hz_h_50m.txt'], [z.' wmode],'delimiter','\t','precision',10);

