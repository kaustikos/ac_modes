close all;
clear variables;
clc;


set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial'); 

cw = 1500;
cb = 1700;

rhow = 1;
rhob = 1.5;

betab = 0.5;

h0 = 200;
zs = 100;
zr = 30;

freq = 25;

omeg = 2*pi*freq;

dz0 = 0.5;
opts.Tgr = 3;
opts.Ngr = 3;
opts.nmod = 8;
opts.Hb = 1000;
opts.BotBC = 'D';

dhs = 180:-2:-166;



kh = [];

MP.HydrologyData = [   [0       cw];
                                 [350    cw]];


for ii = 1:length(dhs)
    
    dh = dhs(ii);
    
    disp(dh);
    
    MP.LayersData = [[0    cw  cw  rhow    rhow       0        0];
                            [h0+dh   cw  cb  rhow    rhob   0   betab]
                             ];
    
    
    if ii == 1
        [krs, wmode, dwmode] = ac_modesr(dz0,MP,freq,opts);
        
        nmod = length(krs);
        kh = zeros(length(dhs),nmod); 
        
       
        
        vgr = zeros(size(kh));
        phizs = zeros(size(kh));
        phizr = zeros(size(kh));
        errs = zeros(size(kh));
        
        
        z = dz0*(0:size(wmode,1)-1);
        
    else
        
        [krs, wmode] = ac_modesr(dz0,MP,freq,opts);
    end;
    
    nmod = length(krs);
    
    wnum_im_part = ModesAttCoeffs(dz0,freq,krs,wmode,MP);
    mgv = ModesGroupVelocities(z,freq,krs,wmode,MP);
    errc = ModesAccuracyCheckPekeris(krs,MP,freq);
    
    
    kh(ii,1:nmod) = krs(1:nmod) + 1i*wnum_im_part;
    vgr(ii,1:nmod) = mgv(1:nmod);
    errs(ii,1:nmod) = errc(1:nmod);
    
    izs = find(z>=zs,1,'first');
    izr = find(z>=zr,1,'first');
    

    
    phizs(ii,1:nmod) = wmode(izs,1:nmod);
    phizr(ii,1:nmod) = wmode(izr,1:nmod);
    
end;


dlmwrite('case_1/kj_wedge_att.txt',[h0+dhs.' kh],'delimiter','\t','precision',10);
dlmwrite('case_1/phizs_wedge.txt',[h0+dhs.' phizs],'delimiter','\t','precision',10);
dlmwrite('case_1/phizr_wedge.txt',[h0+dhs.' phizr],'delimiter','\t','precision',10);
dlmwrite('case_1/vgr.txt',[h0+dhs.' vgr],'delimiter','\t','precision',10);
dlmwrite('case_1/err_pek.txt',[h0+dhs.' errs],'delimiter','\t','precision',10);

figure;
hold all;
for ii = 1:size(kh,2);
    plot(h0+dhs,kh(:,ii));
    
end;


figure;
hold all;
for ii = 1:size(kh,2);
    plot(h0+dhs,phizr(:,ii));
    
end;