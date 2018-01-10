

function out = getFigspecs(figspecs)

[outvars, outvals] = figspechelper(figspecs);
out = cell2struct(outvals',  outvars);

end


function [outvars, outvals] = figspechelper(figspecs);

switch figspecs
    case 'itpc1'
        %ITPC plots
        usecolormap = flipud(cbrewer('div', 'RdBu', 255, 'linear')); %red high white neutral blue low
        contourRes = 1000;
        position = [.1 .1 .9 .8];
        SpacingHorizontal = 0.01;
        SpacingVertical = 0.02;
        Spacing = 0.00;
        Padding = 0.0;
        MarginLeft = 0.04;
        MarginRight = 0.04;
        MarginTop = 0.09;
        MarginBottom =  0.08;
        
    case 'itpc2'
        %ITPC plots
        mcLineColor = [.2 .2 .2];
        usecolormap = 'jet';
        contourRes = 200;
        position = [.1 .1 .9 .8];
        SpacingHorizontal = 0.01;
        SpacingVertical = 0.02;
        Spacing = 0.00;
        Padding = 0.0;
        MarginLeft = 0.04;
        MarginRight = 0.04;
        MarginTop = 0.09;
        MarginBottom =  0.08;
        
    case 'itpcAreas'
        %ITPC plots
        mcLineColor = [.85 .85 .85];
        coloraxis = 'auto';
        usecolormap = 'parula';
        contourRes = 200;
        position = [.1 .1 .7 .3];
        SpacingHorizontal = 0.02;
        SpacingVertical = 0.04;
        Spacing = 0.00;
        Padding = 0.0;
        MarginLeft = 0.05;
        MarginRight = 0.05;
        MarginTop = 0.14;
        MarginBottom =  0.1;
        areatitlesize = 12;
        xticksize = 8;
        
    case 'AreasByDays'
        %ITPC plots
        mcLineColor = [.85 .85 .85];
        coloraxisITPC = [-10 10];
        coloraxisISPC = [-3 3];
        coloraxispower = [-3 3];
        usecolormap = 'parula';
        contourRes = 200;
        position = [0 0 1 1]; %[.1 .1 .7 .8]
        SpacingHorizontal = 0.007;
        SpacingVertical = 0.007;
        Spacing = 0.00;
        Padding = 0.00;
        MarginLeft = 0.05;
        MarginRight = 0.05;
        MarginTop = 0.09;
        MarginBottom =  0.07;
        daytitlesize = 10;
        areatitlesize = 12;
        xtickrotation = 30;
        yticksize = 8;
        xticksize = 8;
        fontname = 'Arial';
        rippleLineWidth = 1;
        daycolor = 'k';
end
% this will put all the variable names in one cell array and their values
% in another. -pats self on back-
iout = whos;
outvars = {iout(:).name};
eval(['outvals = {', strjoin(outvars, ','), '};']);
end