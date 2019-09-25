function [xp, isi]= auxLoadModels(num)

global adaptest
loadVar('.','adaptest',num(1));
m= adaptest{num(1)}{num(2)}{num(3)}{num(4)}.model;
switch m.name
case 'PosPhase_Isi'
    xp= m.xp;
    isi= m.isi;
otherwise
    error(['unknown model ' m.name]);
end
