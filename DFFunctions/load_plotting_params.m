
function plot_params = load_plotting_params(param_set)
%% param_set : cell array of strings
% loads variables for each case of the input param_set to caller workspace
% i.e.:: load_plotting_params({'occnormfiring', 'riptrigspiking'})
% succesors will take priority when var names clash

if ~isa(param_set,'cell')
    param_set = {param_set};
end
for s = param_set
    switch s{1}
        case 'occnormfiring'
            colorSet = 'DR1';
            fonttype = 'Arial';
            titlesize = 16;
            axissize = 16;
            arrowsize = 12;
            usecolormap = 'parula';
            trajname = {'outbound A', 'inbound A', 'outbound B', 'inbound B'};
            SpacingHoriz = 0.02;
            SpacingVert = 0.02;
            Padding = 0.00;
            position = [.1 .1 .8 .5];
            Spacing = 0.00;
            Padding = 0.00;
            MarginLeft = 0.05;
            MarginRight = 0.05;
            MarginTop = 0.15;
            MarginBottom =  0.1;
        case 'riptrigspiking'
            colorSet = 'DR1';
            subtitlecolor = [.6 .6 .6];
            MarkerFaceColor = [.7 .7 .7];
            MarkerFaceAlpha = .7;
            markersize = 10;
            ripMarkerEdgeColor = ([130, 214, 130])./255;
            ripMarkerEdgeAlpha = 1;
            lineColor = ([117, 87, 130])./255;
            lineAlpha = .8;
            siglineColor = ripMarkerEdgeColor;
            siglineAlpha = .5;
            lineWidth = 3;
            areaFaceColor = [.2 .2 .2]; lineColor;
            areaAlpha = lineAlpha;
            SpacingHoriz = 0.02;
            SpacingVert = 0.0;
            Padding = 0.00;
            position = [1.1 1.1 .3 .5];
            Spacing = 0.00;
            Padding = 0.00;
            MarginLeft = 0.05;
            MarginRight = 0.01;
            MarginTop = 0.1;
            MarginBottom =  0.05;
            % create smoothing kernel
            binsize = .001; %% size of bins (in sec)
            % make sure this matches the binsize in filter params
            smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
            smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
            kernel = gaussian(smoothing_width,smoothing_width*8);
    end
end

w = whos;
for a = 1:length(w)
plot_params.(w(a).name) = eval(w(a).name);
end
end