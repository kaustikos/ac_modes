function dscDepths = MPFindDsc(MP)



if isstruct(MP)
    MediaParams = MP.LayersData;
else
    MediaParams = MP;
end;


dscTol = 0.005;

nMP = length(MediaParams(:,1));
ddepths(1:nMP) = MediaParams(1:nMP,1);
cm(1:nMP) = MediaParams(1:nMP,2);
cp(1:nMP) = MediaParams(1:nMP,3);
dm(1:nMP) = MediaParams(1:nMP,4);
dp(1:nMP) = MediaParams(1:nMP,5);

isDsc = ( (abs(cm-cp)>dscTol*cm) | (abs(dm-dp)>dscTol*dm) );
dscDepths = ddepths(isDsc);

