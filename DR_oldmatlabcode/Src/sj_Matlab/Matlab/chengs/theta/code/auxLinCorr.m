function [Cmax, phimax]= auxLinCorr(x,phase)

% function [Cmax, phimax]= auxLinCorr(x,phase)
% shift phases such that |C| becomes maximal
stepsize= 0.1;

validind= find(isfinite(phase) & isfinite(x));
x= x(validind);
phase= phase(validind);
if(length(unique(x))<=1 | length(unique(phase))<=1); 
    Cmax= nan;
    phimax= nan;
    return;
end
Cmax= corrcoef(x,phase);
Cmax= Cmax(1,2);
phimax= 0;
for phi=stepsize:stepsize:2*pi
    C= corrcoef(x,mod(phase+phi,2*pi));
    C= C(1,2);
    if abs(C) > abs(Cmax) | (isnan(Cmax) & ~isnan(C))
        Cmax= C;
        phimax= phi;
    end
end
