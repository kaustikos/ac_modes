function IntRes = CoefIntegrationPiecewise(ziDsc, M1yz, M2zj, dz)

% computation of an integral of the type \int \alpha\phi_i dz,
% where \alpha is a coefficient, which may depend on y and z (M1yz), say \gamma, \nu, etc
% and \phi is some expression consisting of normal modes and their
% derivatives, which depends on z and j (mode number). M1 is given as an
% array ny x nz, M2 as nz x nm. ziDsc contains indices of discontinuities
% in z.


if ~isempty(ziDsc)
    if ziDsc(1)>1
        ziDscNoSurface = ziDsc;
    else
        ziDscNoSurface(1:length(ziDsc)-1) = ziDsc(2:end);
    end;
    
    % IntRes = dz*(  M1yz*M2zj - 0.5*M1yz(:,1)*M2zj(1,:) - 0.5*M1yz(:,end)*M2zj(end,:)+...
    %     0.5*M1yz(:,ziDscNoSurface)*M2zj(ziDscNoSurface,:) - 0.5*M1yz(:,ziDscNoSurface-1)*M2zj(ziDscNoSurface-1,:)   );
    
    IntRes = dz*(  M1yz*M2zj - 0.5*M1yz(:,1)*M2zj(1,:) - 0.5*M1yz(:,end)*M2zj(end,:)-...
        0.5*M1yz(:,ziDscNoSurface)*M2zj(ziDscNoSurface,:) + 0.5*M1yz(:,ziDscNoSurface-1)*M2zj(ziDscNoSurface-1,:)   );
else
    IntRes = dz*(  M1yz*M2zj - 0.5*M1yz(:,1)*M2zj(1,:) - 0.5*M1yz(:,end)*M2zj(end,:));
end;