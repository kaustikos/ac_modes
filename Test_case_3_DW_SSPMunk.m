%close all;
clear variables;
clc;



MP.LayersData = [[0 1500 1500 1 1 0 0] ];

dz = 0.5;
z = 0:dz:3000;
freq = 400;
zbar = 2*(z-1300)/1300;
eps = 0.00737;
c = 1500*(1 + eps*(zbar - 1 + exp(-zbar)));


MP.HydrologyData = [z.' c.'];


opts.nmod = 30;
opts.Hb = 3000;
opts.Ngr = 3;
opts.BotBC = 'D';
opts.Tgr = 3;



tic
%[wnum, wmode] = ac_modes(z,MP,freq,0.7,'D');

[krs, wmode] = ac_modesr(dz,MP,freq, opts );

disp(krs.');
toc


z = dz*(0:size(wmode,1)-1);

figure;
hold all;
for ii = 1:5:opts.nmod
    
    [~,im] = max(abs(wmode(:,ii)));
    
    wmode(:,ii) = wmode(im,ii)*wmode(:,ii);
    
    plot(z,  wmode(:,ii));
end;



mgv = ModesGroupVelocities(z,freq,krs,wmode,MP);

figure;
plot(krs+1i*0.00000001,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);


dlmwrite(['case_3/kj_f_' int2str(freq) 'Hz.txt'], (krs).','delimiter','\t','precision',10);
dlmwrite(['case_3/vgj_f_' int2str(freq) 'Hz.txt'], (mgv).','delimiter','\t','precision',10);
dlmwrite(['case_3/phij_f_' int2str(freq) 'Hz.txt'], [z.' wmode],'delimiter','\t','precision',10);

