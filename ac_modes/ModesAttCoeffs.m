function wnum_im_part = ModesAttCoeffs(dz,freq,wnum,wmode,MP)

omeg = 2*pi*freq;


nz = size(wmode,1);

z = (0:nz-1)*dz;

[ziDsc, c, d, beta] = MediaParamsToVectors(z,MP);

eta = 1/(40*pi*log10(exp(1)));

wnum_im_part = (omeg^2*eta)*CoefIntegrationPiecewise(ziDsc, beta./( d.*(c.^2) ), wmode.^2, dz)./wnum;


