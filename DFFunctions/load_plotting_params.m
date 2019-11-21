
function Pp = load_plotting_params(param_set, varargin)
% Pp = load_plotting_params(param_set, varargin)
% param_set : cell array of strings
% each string matches a param set
% succesors take priority if vars clash
%{
FFPhy V0.1
@DR
%}

if ~isempty(varargin)
    assign(varargin{:})
end
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
    MgRt = 0.02; % margin right
    MgTp = 0.08; % margin top
    MgBm =  0.08; % margin bottom
    FontS = 8;
    FontW = 'bold';
    FontNm = 'Arial';
    SupFontS = 12;
    tickSz = 8; % tick labels font size
    set(0,'defaultAxesFontSize',10)


%% ========= SWR 'eventTrigSpiking'=========

case 'wtrackSWRSU'
    position = [.1 .1 .25 .4];
    bin = .010; %sec
    win = [.5 .5];
    SpVt = 0.01;
    MgLt = 0.1; % margin left
    MgBm = 0.1; % margin bottom
    spikeSz = 50; 
    spikeAlpha = .5;
    
case 'wtrackSWRHeatrast'
    position = [.1 .1 .8 .5];
    bin = .010; %sec
    win = [.3 .3];
    SpHz = 0.05;
    MgLt = 0.1; % margin left
    MgBm = 0.1; % margin bottom
    
case 'wtrackSWRHeatRastAllAn'
    position = [.1 .1 .3 .5];
    bin = .010; %sec
    win = [.3 .3];
    SpHz = 0.05;
    MgLt = 0.1; % margin left
    MgBm = 0.1; % margin bottom

case 'wtrackSWRModCDF'
    
%% ========= Lick 'eventTrigSpiking'=========

case 'wtrackLickSU'
    position = [.1 .1 .25 .4];
    bin = .010; %sec
    win = [.5 .5];
    SpVt = 0.01;
    MgLt = 0.1; % margin left
    MgBm = 0.1; % margin bottom
    spikeSz = 50; 
    spikeAlpha = .5;
    
case 'wtrackLickHeatrast'
    position = [.1 .1 .3 .5];
    bin = .010; %sec
    win = [.3 .3];
    SpHz = 0.05;
    MgLt = 0.1; % margin left
    MgBm = 0.1; % margin bottom    
    
case 'wtrackLickHeatRastAllAn'
    position = [.1 .1 .3 .5];
    bin = .010; %sec
    win = [.3 .3];
    SpHz = 0.05;
    MgLt = 0.1; % margin left
    MgBm = 0.1; % margin bottom

case 'wtrackLickModCDF'

%% ========= 
    
case 'SWRTrigHeatraster'
    position = [.1 .1 .3 .5];
    bin = .010; %sec
    win = [.5 .5];
    SpHz = 0.05;
    MgLt = 0.1; % margin left
    MgBm = 0.1; % margin bottom

case 'SULickILPC_cdfPolar'
    position = [.1 .1 .2 .4];
    logZBins = 200;
    sigPct = .95;
    SpVt = 0.15;
    MgLt = 0.15;

case 'suILPC'
    position = [.1 .1 .2 .4];
    MgLt = 0.15;
    SpVt = 0.01;
    numBins = 16;
    bin = .010; %sec
    win = [.5 .5];
    SpHz = 0.05;

case 'spikeILPC'
    position = [.1 .1 .5 .7];
    bin = .010; %sec
    win = [.5 .5];
    SpHz = 0.05;

case 'ILIvIBI'
    position = [.1 .1 .5 .6];
    SpHz = 0.1; %spacing horizontal
    SpVt = 0.1; % spacing vertical

case 'HeatRastSUmod'
    position = [.1 .1 .5 .7];
    bin = .010; %sec
    win = [.5 .5];
    SpHz = 0.05;

case 'lickTrigSUmod'
    position = [.1 .1 .4 .7];
    spikeSz = 50;
    win = [.5 .5];
    spikeAlpha = .6;
    bin = .010; %sec
    SpVt= 0;

case 'lickTrigSpikingPSTH'
    position = [.1 .1 .5 .5];
    psthSize = 5; 

case 'RxnFull'
    position = [.1 .1 .5 .5];
    winSE = [-.5 .5];
    interpBy = 10;

case 'RxnPerPC'
    position = [.1 .1 .5 .5];
    winSE = [-.5 .5];

case 'areaTFspect'
%             area = {{'mec', 'sup'}, {'mec', 'deep'}, {'ca1', 'd'}};
    position = [.1 .1 .8 .5];
    win = [-.5 .5];
    usecolormap = hot;
    tickFsize = 8;

case 'wetVDryILIphaseSWR'
    position = [.1 .1 .5 .5];

case 'FeatureTracking'
    SpVt = 0.00;
    Msz = 30;
    lickClr = [.5 0 1 .2];
    rewClr = [.3 .3 .3];
    position = [.1 .1 .3 .4];

case 'lickBoutSUswrmod'
    position = [.1 .1 .6 .7];
    SpVt = 0.05;
    win = [.5 .5];
case 'examples'
    position = [.1 .1 .6 .5];
    Yoffset = 600;
    nrows = 5;
    SpHz = 0;
    MgBm =  0.1;
    MgRt =  0.04;
    MgLt = .08;
    MgTp = .1;
case 'fitLM'
    stitFsize = 16;
    tickFsize = 8;
    sfTitFsize = 14;
    pwin = [1 1];
    usecolormap = 'jet';
    position = [.1 .1 1 1];
    SpHz = 0.03; %spacing horizontal
    SpVt = 0.05; % spacing vertical
    MgLt = 0.04; % margin left
    MgTp = 0.07; % margin top
    MgBm =  0.04; % margin bottom

case 'dfa_lickswrcorr'
    position = [.1 .1 .5 .5];
    SpHz = 0.15; %spacing horizontal
    SpVt = 0.15; % spacing vertical

case 'xcorrPhasemodSU'
    sigpct = 97.5; % pct
    position = [.1 .1 .6 .6];
    MgBm =  0.07; % margin bottom
    SpVt = 0.1;
    SpHz = 0.1;
    MgLt = .11;
    MgTp = .06;
    psthSize = 5;
    win = [.5 .5];

case 'dfa_lickXCorrSpikes'
    position = [.1 .1 .3 .5];
    MgBm =  0.07; % margin bottom
    SpVt = 0.1;
    SpHz = 0.1;
    MgLt = .11;
    MgTp = .06;
    psthSize = 5;
    win = [.5 .5];

case 'combinedAreasTFstats'
    stitFsize = 16;
    tickFsize = 8;
    sfTitFsize = 14;
    pwin = [1 1];
    usecolormap = 'bone';
    contourRes = 40;
    position = [.1 .1 .5 .7];
    SpHz = 0.04; %spacing horizontal
    SpVt = 0.06; % spacing vertical
    MgLt = 0.02; % margin left
    MgTp = 0.1; % margin top
    MgBm =  0.06; % margin bottom
    MgLt = 0.1; % margin left
case 'powerTFmap'
    stitFsize = 16;
    tickFsize = 8;
    sfTitFsize = 14;
    pwin = [.5 .5];
    usecolormap = 'bone';
    contourRes = 40;
    position = [.1 .1 1 1];
    SpHz = 0.02; %spacing horizontal
    SpVt = 0.04; % spacing vertical
    MgLt = 0.02; % margin left
    MgTp = 0.1; % margin top
    MgBm =  0.04; % margin bottom
%         case 'expvarCatMeanPwr'
% %             mcLineColor = [.2 .2 .2];
%             pwin = [1 1];
%             usecolormap = 'jet';
%             contourRes = 40;
%             position = [.1 .1 1 1];
%             SpHz = 0.02; %spacing horizontal
%             SpVt = 0.04; % spacing vertical
%             MgLt = 0.02; % margin left
%             MgTp = 0.1; % margin top
%             MgBm =  0.04; % margin bottom


case 'behaveperform'
    colorSet = 'DR1';
    % clims = [0 1]; %[0 .7]
    % position = [.1 .1 .9 .8];
    SpacingHorizontal = 0.01;
    SpacingVertical = 0.02;
    % Spacing = 0.00;
    % Padding = 0.0;
    % MarginLeft = 0.04;
    % MarginRight = 0.04;
    % MarginTop = 0.09;
    % MarginBottom =  0.08;
    position = [.1 .1 .5 .7];
    SpacingHorizontal = 0.00;
    SpacingVertical = 0.00;
    Spacing = 0.00;
    Padding = 0.00;
    MarginLeft = 0.05;
    MarginRight = 0.05;
    MarginTop = 0.14;
    MarginBottom =  0.08;
    usecolormap = 'jet';
    win = [.5 .5]; %in seconds
    % indwin = win*1500;
case 'frmaps'
    position = [.1 .1 .8 .4];
    tits = {'outbound right', 'inbound right', 'outbound left', 'inbound left'};
%             Sp = 0.00; % spacing
    SpHz = 0.02; %spacing horizontal
    SpVt = 0.02; % spacing vertical
%             Pad = 0.00; % padding
    MgLt = 0.05; % margin left
    MgRt = 0.05; % margin right
    MgTp = 0.1; % margin top
    MgBm =  0.1; % margin bottom
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
    position = [.1 .1 .9 .3];
    pwin = [1 1];
case 'powerheatRast'
    position = [.1 .1 .9 .9];
    pwin = [1 1];
    plot_frex = [8 16 40 85 150 300];
    cmap = 'bone';
case 'riptrigvel'
    position = [.1 .1 1 .5];
    SpHz = 0.02;
    SpVt = 0.08;
    MgBm = .1;
case 'dfa_plotDataChunks'
    position = [0 0 1 1];
    nrows = 15;
    ripFAlpha = .4;
    threshFAlpha = .5;
    SUsize = 6;
    areaSepWidth = 2;
    areaSepAlpha = 1;
    rippwrcmap = 'hot';
    thetapwrcmap = 'parula';
    ripclr = 'y';
    noboutsclr = 'b';
    boutsclr = 'r';
    noboutsalpha = .05;
    boutsalpha = .05;

    c1clr = colorPicker({'ca1'},'subtags', {'d'});
    dmclr = colorPicker({'mec'},'subtags', {'deep'});
    smclr = colorPicker({'mec'},'subtags', {'supf'});
    phaseCmap = 'pink';
    inputclr = colorPicker({'well'}, 'subtags', {'input'});
    outputclr = colorPicker({'well'}, 'subtags', {'output'});
    dioEAlpha = 0;
    dioFAlpha = .25;
    linposclr = [.3 .3 .3];
    linposcmap = cbrewer('qual', 'Set2', 5);
    linposSegAlpha = .25;
    bpzrippwr = [0 .8 1 .7];
%             SpHz = 0.005; %spacing horizontal
    SpVt = 0; % spacing vertical
    MgTp = 0.04;
    MgLt = 0.06;
    MgRt = 0.0;
    MgBm = 0.05;
case 'riptriglfp_allLFPtype_perntrode'
    position = [.1 .1 .4 .7];
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
%         case 'statespacePerformance'
%             position = [.05 .05 .5 .5];
case 'riptrigFRxdays'
    position = [.05 .05 .3 .6];
end
end

w = whos;
for a = 1:length(w)
Pp.(w(a).name) = eval(w(a).name);
Pp.posparams = {'SpacingVert', SpVt, 'SpacingHoriz', SpHz, ...
    'MarginLeft', MgLt, 'MarginRight', MgRt, 'MarginTop', MgTp, ...
    'MarginBottom', MgBm};
end
% if init_fig
%     init_plot(showfigs, 'position', plot_params.position)
% end
end