
function load_plotting_params(param_set)
%% param_set : cell array of strings
% loads variables for each case of the input param_set to caller workspace
% i.e.:: load_plotting_params({'occnormfiring', 'riptrighist'})

if ~isa(param_set,'cell')
    param_set = {param_set};
end
for s = param_set
    switch s{1}
        case 'occnormfiring'
            P.colorSet = 'DR1';
            P.fonttype = 'Arial';
            P.titlesize = 16;
            P.axissize = 16;edit 
            P.arrowsize = 12;
            P.usecolormap = 'hot';
            P.trajname = {'outbound A', 'inbound A', 'outbound B', 'inbound B'};
            P.SpacingHoriz = 0.02;
            P.SpacingVert = 0.02;
            P.Padding = 0.00;
            P.position = [.1 .1 .8 .5];
            P.Spacing = 0.00;
            P.Padding = 0.00;
            P.MarginLeft = 0.05;
            P.MarginRight = 0.05;
            P.MarginTop = 0.15;
            P.MarginBottom =  0.1;
    end
end
   assignin('caller', 'P', P)
end