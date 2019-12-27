function [wnum, wmode, varargout] = ac_modesr(dz0,MediaParams,freq, varargin)

% Solution of acoustical spectral problem + Richardson extrapolation;

%% options

% number of grids to solve the problem

Ngr = 3;

% grid type controls mesh sizes for extrapolation
% mesh size factors: 1 = half-integers (1, 1.5, 2); 2 = integers (1,2,3); 
% 3 = fixed sizes, 0.25 m, 0.5 m, 1 m

Tgr = 1;

% fictitious bottom depth

if isstruct(MediaParams)
    LayersData = MediaParams.LayersData;
else
    LayersData = MediaParams;
end

if LayersData(end,1) >= 500
    Hb = LayersData(end,1) + 500; % deep water: +500 m below the last LayersData entry <-heuristic!
elseif LayersData(end,1) < 100
    Hb = 200; % shallow water: fictitious bottom at 100 m;
else
    Hb = 2*LayersData(end,1);
end

% nmod: trapped modes only = -1 (default); propagating modes = 0; given number > 0

nmod = -1;
BotBC = 'D';

if nargin == 4
    if isstruct(varargin{1})
        opts = varargin{1};
        
        if isfield(opts,'nmod')
            nmod = opts.nmod;
        end
        if isfield(opts,'Ngr')
            Ngr = opts.Ngr;
        end
        if isfield(opts,'Tgr')
            Tgr = opts.Tgr;
        end
        if isfield(opts,'Hb')
            Hb = opts.Hb;
        end
        if isfield(opts,'BotBC')
            BotBC = opts.BotBC;
        end
    else
        nmod = varargin{1};
    end
end

% find bottom depth (the first discontinuity)
% if it doesn't exit than use fictitious bottom depth

dscDepths = MPFindDsc(MediaParams);

if ~isempty(dscDepths)
    bDepth = dscDepths(1);
else
    bDepth = Hb;
end



%% wavenumbers computation: loop over Richardson grids

% Solving spectral problem for several grids with different meshsize
% initializing arrays 

wnums(1:Ngr,1:20) = 0; % wnumbers for different meshsizes
dzs(1:Ngr) = dz0; % meshsizes

for ii = 1:Ngr
     
    if Tgr == 1
        izH = fix( bDepth/(dz0*(0.5+ii/2) ) ); 
    elseif Tgr == 2
        izH = fix( bDepth/(dz0*(ii) ) );    
    end
    
    
    if Tgr < 3
        if izH > 0
            dzc = bDepth/izH;
        else
            error(['ac_modesr: Water layer is too thin for a given grid. Ngr=' num2str(Ngr) '; bdepth=' num2str(bDepth)]); 
        end
    else
        dzc = dz0*(2^(ii-1));
    end
    zc = 0:dzc:Hb;
    
    dzs(ii) = dzc;    
    
    if ii==1
        if nargout < 3
            [wnum, wmode] = ac_modes(zc,MediaParams,freq, nmod, BotBC);
        else
            [wnum, wmode, dwmode] = ac_modes(zc,MediaParams,freq, nmod, BotBC);
            varargout{1} = dwmode;
        end
        nmodc = length(wnum);    
    else
        [wnum, ~] = ac_modes(zc,MediaParams,freq, nmod, BotBC);
        
        nmodc = min(nmodc, length(wnum));
        
    end
    
    wnums(ii,1:nmodc) = wnum(1:nmodc);    
end

% Richardson extrapolation

rMatrix(1:Ngr,1:Ngr) = repmat(dzs(1:Ngr).',1,Ngr).^repmat( 2*(0:Ngr-1),Ngr,1);
rDat = rMatrix\(wnums(1:Ngr,1:nmodc).^2);
wnum = sqrt( rDat(1,1:nmodc) );

%AA = rMatrix^(-1);
%disp('Richardson coeffs:');
%disp( AA(1,:) );


