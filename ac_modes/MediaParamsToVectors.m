function [ziDsc, varargout] = MediaParamsToVectors(z,MediaParams)

%syntax1: ziDsc = MediaParamsToVectors(z,MediaParams) - construct array of
%z-indices of media discontinuities

%syntax2: [ziDsc, c, d] = MediaParamsToVectors(z,MediaParams) - construct array of
%z-indices of media discontinuities, sound speed and density vectors
%(possibly having finite jumps at z in ziDsc)

%syntax3: [ziDsc, c, d, h0] = MediaParamsToVectors(z,MediaParams) - construct array of
%z-indices of media discontinuities, sound speed and density vectors
%(possibly having finite jumps at z in ziDsc) and interface depths h0

nz = length(z);
dz = z(2) - z(1);

h0 = MPFindDsc(MediaParams);
ziDsc = round(h0/dz)+1;


if isstruct(MediaParams)
    LayersData = MediaParams.LayersData;
else
    LayersData = MediaParams;
end;


nDsc = length(LayersData(:,1));
ziParams(1:nDsc) = round(LayersData(1:nDsc,1)/dz)+1;
ddepths(1:nDsc) = LayersData(1:nDsc,1);
cm(1:nDsc) = LayersData(1:nDsc,2);
cp(1:nDsc) = LayersData(1:nDsc,3);
dm(1:nDsc) = LayersData(1:nDsc,4);
dp(1:nDsc) = LayersData(1:nDsc,5);

%% constructing sound speed profile

if nargout>1
c(1:nz)=cm(1);

for ii=2:nDsc
    if ziParams(ii-1)~=ziParams(ii)
        c(ziParams(ii-1):ziParams(ii))=interp1([ddepths(ii-1), ddepths(ii)],...
            [cp(ii-1), cm(ii)],z(ziParams(ii-1):ziParams(ii)),'linear','extrap');
    else
        c(ziParams(ii)) = cp(ii);
    end;
end;

c(ziParams(nDsc):end)=LayersData(nDsc,3);

if isstruct(MediaParams)
    if isfield(MediaParams,'HydrologyData')
        if nDsc > 1
            if ziParams(2) > 1
                c(1:ziParams(2)) = interp1(MediaParams.HydrologyData(:,1),MediaParams.HydrologyData(:,2),z(1:ziParams(2)),'linear','extrap');
                %???cm(2) = c(ziParams(2));
                c(ziParams(2)) = cp(2);
            end;
        else
            c = interp1(MediaParams.HydrologyData(:,1),MediaParams.HydrologyData(:,2),z,'linear','extrap');
        end;
    end;
end;

    varargout{1} = c;
end;


%% constructing density profile

if nargout>2
    
    d(1:nz)=dm(1);
    
    for ii=2:nDsc
        if ziParams(ii-1)~=ziParams(ii)
            d(ziParams(ii-1):ziParams(ii))=interp1([ddepths(ii-1), ddepths(ii)],...
                [dp(ii-1), dm(ii)],z(ziParams(ii-1):ziParams(ii)),'linear','extrap');
        else
            d(ziParams(ii)) = dp(ii);
        end;
    end;
    
    d(ziParams(nDsc):end)=LayersData(nDsc,5);
   
    varargout{2} = d;
end;

%% constructing attenuation profile

if nargout>3
    attm(1:nDsc) = LayersData(1:nDsc,6);
    attp(1:nDsc) = LayersData(1:nDsc,7);
    
    att(1:nz)=attm(1);
    
    for ii=2:nDsc
        if ziParams(ii-1)~=ziParams(ii)
            att(ziParams(ii-1):ziParams(ii))=interp1([ddepths(ii-1), ddepths(ii)],...
                [attp(ii-1), attm(ii)],z(ziParams(ii-1):ziParams(ii)),'linear','extrap');
        else
            att(ziParams(ii)) = attp(ii);
        end;
        
    end;
    

    att(ziParams(nDsc):end)=LayersData(nDsc,7);
    varargout{3} = att;
end;

%% depths

if nargout>4
    varargout{4} = h0;
end;