function mgv = ModesGroupVelocities(z,freq,wnum,wmode,MP)

% warning! this code doesn't work well! some error is here
% 11.08.2017

omeg = 2*pi*freq;

dz = z(2) - z(1);
nM = length(wnum);
nz = length(z);

[ziDsc, c, d] = MediaParamsToVectors(z,MP);
imgv = omeg*CoefIntegrationPiecewise(ziDsc, 1./( d.*(c.^2) ), wmode.^2, dz)./wnum;

mgv = 1./imgv;

