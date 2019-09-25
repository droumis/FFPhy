function var= getCellSelectVar(vname)

if nargin < 1
    error('no variable name given');
end

global fmaux select

if isfield(select, vname)
    var= select.(vname)(fmaux.currentCell);
else
    warning([vname ' is not a variable in select ' fmaux.select...
	     '. Setting to default.']);
    if strcmp(vname, 'traj')
        var= [0:3];
    end
end

