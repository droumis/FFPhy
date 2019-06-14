
function plot_params = load_plotting_params(param_set)
%% param_set : cell array of strings
% loads variables for each case of the input param_set to caller workspace
% i.e.:: load_plotting_params({'occnormfiring', 'riptrigspiking'})
% succesors will take priority when var names clash

% helpful: 
% 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', Pp.FontNm

if ~isa(param_set,'cell')
    param_set = {param_set};
end
for s = param_set
    switch s{1}
        case 'defaults'
            position = [.1 .1 .4 .4];
            pwin = [1 1];
            %             Sp = 0.00; % spacing
            SpHz = 0.005; %spacing horizontal
            SpVt = 0.05; % spacing vertical
            %             Pad = 0.00; % padding
            MgLt = 0.08; % margin left
            MgRt = 0.01; % margin right
            MgTp = 0.1; % margin top
            MgBm =  0.05; % margin bottom
            FontS = 8;
            FontW = 'bold';
            FontNm = 'Arial';
            SupFontS = 12;
            tickSz = 8; % tick labels font size
        case 'riptriglfp_perntrode_perep_traces'
            position = [.1 .1 .4 .4];
            pwin = [1 1];
%             Sp = 0.00; % spacing
            SpHz = 0.06; %spacing horizontal
            SpVt = 0.02; % spacing vertical
%             Pad = 0.00; % padding
            MgLt = 0.1; % margin left
            MgRt = 0.05; % margin right
            MgTp = 0.1; % margin top
            MgBm =  0.1; % margin bottom
            FontS = 8;
            FontW = 'bold';
            FontNm = 'Arial';
            SupFontS = 12;
            tickSz = 8; % tick labels font size
        case 'ripcleaning'
            position = [.05 .05 .75 .75];
            pre_excl_win = 30; %seconds
            post_excl_win = 2; %seconds
            ax_srate = 30; %Hz
            cm_thresh = 10; %cm
            sig = 400;
            patchtxtsz = 8;
        case 'riptrigall'
            position = [.05 .05 .75 .75];
            patchtxtsz = 8;
            stdthresh = 2;
            t = 10;
            tscale = [1000 1000 t t t t t];
            nrow = 8;
            SpacingHoriz = 0.02;
            SpacingVert = 0.04;
            Padding = 0.00;
            Spacing = 0.00;
            Padding = 0.00;
            MarginLeft = 0.05;
            MarginRight = 0.05;
            MarginTop = 0.1;
            MarginBottom =  0.1;
        case 'riptriglfp_perLFPtype_allntrodes'
            position = [.1 .1 .9 .9];
            pwin = [1 1];
        case 'riptriglfp_allLFPtype_perntrode'
            position = [.1 .1 .9 .4];
            pwin = [1 1];
        case 'riptriglfp_acrossdays'
            position = [.1 .1 .8 .3];
        case 'riptriglfp_acrossdays_allntrodes'
            position = [.1 .1 .8 .8];
            
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
            position = [1 1 1 .4];
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
        case 'all_nts_days_riptrigspiking'
            position = [.05 .05 .95 .95];
            
        case 'ripcorrVperform_wriptrig'
            position = [.05 .05 .95 .95];
        case 'ripcorr_Xdays'
            position = [.05 .05 .5 .5];
            cmap = colormap(lines);
        case 'paircorrVperformance'
            position = [.05 .05 .95 .95];
            alpha = .001;
        case 'paircorrVperformance_hm'
            position = [.05 .05 .5 .5];
        case 'statespacePerformance'
            position = [.05 .05 .5 .5];
        case 'riptrigFRxdays'
            position = [.05 .05 .3 .6];
    end
end

w = whos;
for a = 1:length(w)
plot_params.(w(a).name) = eval(w(a).name);
end
end