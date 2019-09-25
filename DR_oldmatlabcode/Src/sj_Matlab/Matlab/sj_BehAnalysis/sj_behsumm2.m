
function [summ1] = sj_behsumm2(prefixes, expidx, conidx, normidx, figcurve, figsmooth, figout, figin, savefig1,savedata1)
% sj_behsumm2({'REc'; 'REd'; 'REe'; 'REf'; 'REg';'REh';'RCa';'RCb';'RCc';'RCd';'RNa';'RNb';'RNc';'RNd'},[1:6],[7:10],[11:14],1,0,0,0,0,0);

% sj_behsumm2('/data25/sjadhav/RippleInterruption/Figures/01May11_RippleDisFigs/Behavior',{'REc'; 'REd'; 'REe'; 'REf';'RCa';'RCb';'RCc';'RCd';'RNa';'M25';'M26';'M06'},[1:4],[5:8],[9:12],1,0,0)
% sj_behsumm2('/data25/sjadhav/RippleInterruption/Figures/01May11_RippleDisFigs/Behavior',{'REc'; 'REd'; 'REe'; 'REf';'RCa';'RCb';'RCc';'RCd';'Cor';'Fiv';'Sev';'Eig';'ten'},[1:4],[5:8],[9:13],1,0,0);
% sj_behsumm2('/data25/sjadhav/RippleInterruption/Figures/01May11_RippleDisFigs/Behavior',{'REc'; 'REd'; 'REe'; 'REf';'RCa';'RCb';'RCc';'RCd'},[1:4],[5:8],[],1,0,0);
% sj_behsumm2('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'REc'; 'REd'; 'REe'; 'REf';'RCa';'RCb';'RCc';'RNa'},[1:4],[5:8],[],1,0,0)
% sj_behsumm2('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'REc'; 'REd'; 'REe'; 'REf';'RCa';'RCb';'M25';'M26'},[1:4],[5:8],[],1,0,0)
% sj_behsumm2('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'REc'; 'REd'; 'REe'; 'REf'; 'RE1';'RCa';'M25';'M26'},[1:4],[5:8],[],1,0,0);
% sj_behsumm2 ('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'RE1';'RCa'; 'REc'; 'REd'; 'REe'; 'REf'; 'dud'},[3:6],[1:2],[7],0,0,0);
% sj_behsumm2('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'REc'; 'REd'; 'REe';'M25';'M26';'RCa';'M06'},[1:3],[4:7],[],1,0,0);

% f=1; % Choice For Plotting All Performance Curves
% f2=0; % Choice For Plotting ot Omitting Smooth Performance Curves
% f3=0; % Choice for plotting Group Comparison Stats

% Make behavior summary for all groups. Smooth_perf, adaptive algorithm, Ntrajs, Speed, etc.
% Shantanu Jadhav, 02/22/11; Post Rips-talk

%

if nargin<5,
    figcurve = 0; % Choice for plotting Fitted Curves
end

if nargin<6,
    figsmooth = 0;  % Choice For Plotting Moving Average Curves, and Raw Perf vs Day
end

if nargin<7,
    figout = 0; % Choice for plotting all Outbound Summary Group Comparison Stats
end

if nargin<8,
    figin = 0; %  Choice for plotting all Inbound Summary Group Comparison Stats
end

if nargin<10,
    savefig1 = 0; % Save figures or not
end

if nargin<11,
    savedata1 = 0; % Save Data or not
end

savefig2=1; % Use manually to save stuff you want

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Behavior/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

clr = {'b','r','g','c','m','y','k','r'};

% ---------------------------------------


% Params
smoothperf_bin=10;
Ntraj_cutoff=200;
day_cutoff=8;
lastdays=2;

% Init
summ1=1;
summdirectoryname = summdir;
if (summdirectoryname(end) == '/')
    summdirectoryname = summdirectoryname(1:end-1);
end
%cd(summdirectoryname);

% Get each Animal

nexp=length(expidx); ncon=length(conidx); nnor=length(normidx);

for n=1:length(prefixes)
    
    currprefix = prefixes{n};
    switch currprefix
        case 'RE1'
            directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct';
        case 'RNa'
            directoryname = '/data25/sjadhav/RippleInterruption/RNa_direct';
        case 'RNb'
            directoryname = '/data25/sjadhav/RippleInterruption/RNb_direct';
        case 'RNc'
            directoryname = '/data25/sjadhav/RippleInterruption/RNc_direct';
        case 'RNd'
            directoryname = '/data25/sjadhav/RippleInterruption/RNd_direct';
        case 'RCa'
            directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct';
        case 'RCb'
            directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct';
        case 'RCc'
            directoryname = '/data25/sjadhav/RippleInterruption/RCc_direct';
        case 'RCd'
            directoryname = '/data25/sjadhav/RippleInterruption/RCd_direct';
        case 'REc'
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct';
        case 'REd'
            directoryname = '/data25/sjadhav/RippleInterruption/REd_direct';
        case 'REe'
            directoryname = '/data25/sjadhav/RippleInterruption/REe_direct';
        case 'REf'
            directoryname = '/data25/sjadhav/RippleInterruption/REf_direct';
        case 'REg'
            directoryname = '/data25/sjadhav/RippleInterruption/REg_direct';
        case 'REh'
            directoryname = '/data25/sjadhav/RippleInterruption/REh_direct';    
        case 'M25'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';
        case 'M26'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';
        case 'M24'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';
        case 'M06'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';
        case 'Cor'
            directoryname = '/data25/sjadhav/RippleInterruption/OthControlData/Cor';
        case 'Fiv'
            directoryname = '/data25/sjadhav/RippleInterruption/OthControlData/Fiv';
        case 'Sev'
            directoryname = '/data25/sjadhav/RippleInterruption/OthControlData/Sev';
        case 'dud'
            directoryname = '/data25/sjadhav/RippleInterruption/OthControlData/Dud';
        case 'Eig'
            directoryname = '/data25/sjadhav/RippleInterruption/OthControlData/Eig';
        case 'ten'
            directoryname = '/data25/sjadhav/RippleInterruption/OthControlData/Ten';
    end
    
    cd(directoryname);
    
    % 1) Outbound behavior
    filename = [currprefix '_outbound'];
    load(filename);
    
    Ripdis_summ(n).outbound_logic=all_outbound_logic;
    Ripdis_summ(n).out_curve=outbound_curve;
    Ripdis_summ(n).out_curve_witherr=outbound_curve_witherr;
    Ripdis_summ(n).out_curve_lowerr=outbound_lowerr;
    Ripdis_summ(n).out_curve_uperr=outbound_uperr;
    Ripdis_summ(n).out_ntrajs_perday=ntrajs_perday_out;
    Ripdis_summ(n).out_ntrajs_persess=ntrajs_persess_out;
    Ripdis_summ(n).out_ntraj_axis_day=ntraj_axis_day;
    Ripdis_summ(n).out_ntraj_axis_sess=ntraj_axis_sess;
    Ripdis_summ(n).out_learning_trial=outbound_learningtrial;
    Ripdis_summ(n).out_learning_sess=learning_sess;
    Ripdis_summ(n).out_learning_day=learning_day;
    Ripdis_summ(n).outbound_perf_day=outbound_perf_day;
    
    [perf, smooth_perf, rawperfvsday] = sj_smoothperf1 (currprefix, smoothperf_bin, all_outbound_logic, ntrajs_perday_out, 'Out',0);
    % Note: perf = smooth_perf. No difference
    
    Ripdis_summ(n).out_smoothperf=smooth_perf;
    Ripdis_summ(n).out_runavg=perf;
    Ripdis_summ(n).out_rawperfvsday=rawperfvsday;
    
    if n<=nexp+ncon
        if ~( strcmp(currprefix,'RCa') | strcmp(currprefix,'REc') )
            Ripdis_summ(n).ex_outbound_perf_day=ex_outbound_perf_day;
            Ripdis_summ(n).ex_out_curve=ex_outbound_curve;
            Ripdis_summ(n).ex_outbound_logic=ex_outbound_logic;
            
            Ripdis_summ(n).expre_outbound_perfvsday=outbound_perf_day;
            Ripdis_summ(n).expre_outbound_rawperfvsday=rawperfvsday;  % same as out_rawperfvsday
            
        else
            Ripdis_summ(n).ex_outbound_perf_day=[];
            Ripdis_summ(n).ex_out_curve=[];
            Ripdis_summ(n).ex_outbound_logic=[];
            
            Ripdis_summ(n).expre_outbound_perfvsday=[];
            Ripdis_summ(n).expre_outbound_rawperfvsday=[];
        end
    end
    
    
    
    % 2) Inbound Behavior
    filename = [currprefix '_inbound'];
    load(filename);
    
    Ripdis_summ(n).inbound_logic=all_inbound_logic;
    Ripdis_summ(n).in_curve=inbound_curve;
    Ripdis_summ(n).in_curve_witherr=inbound_curve_witherr;
    Ripdis_summ(n).in_curve_lowerr=inbound_lowerr;
    Ripdis_summ(n).in_curve_uperr=inbound_uperr;
    Ripdis_summ(n).in_ntrajs_perday=ntrajs_perday_in;
    Ripdis_summ(n).in_ntrajs_persess=ntrajs_persess_in;
    Ripdis_summ(n).in_ntraj_axis_day=ntraj_axis_day;
    Ripdis_summ(n).in_ntraj_axis_sess=ntraj_axis_sess;
    Ripdis_summ(n).in_learning_trial=inbound_learningtrial;
    Ripdis_summ(n).in_learning_sess=learning_sess;
    Ripdis_summ(n).in_learning_day=learning_day;
    Ripdis_summ(n).inbound_perf_day=inbound_perf_day;
    
    [perf, smooth_perf, rawperfvsday] = sj_smoothperf1 (currprefix, smoothperf_bin, all_inbound_logic, ntrajs_perday_in, 'In',0);
    % Note: perf = smooth_perf. No difference
    
    Ripdis_summ(n).in_smoothperf=smooth_perf;
    Ripdis_summ(n).in_runavg=perf;
    Ripdis_summ(n).in_rawperfvsday=rawperfvsday;
    
    if n<=nexp+ncon
        if ~( strcmp(currprefix,'RCa') | strcmp(currprefix,'REc') )  % No extra days for these animals
            Ripdis_summ(n).ex_inbound_perf_day=ex_inbound_perf_day;
            Ripdis_summ(n).ex_in_curve=ex_inbound_curve;
            Ripdis_summ(n).ex_inbound_logic=ex_inbound_logic;
            
            Ripdis_summ(n).expre_inbound_perfvsday=inbound_perf_day;
            Ripdis_summ(n).expre_inbound_rawperfvsday=rawperfvsday;
            
        else
            Ripdis_summ(n).ex_inbound_perf_day=[];
            Ripdis_summ(n).ex_in_curve=[];
            Ripdis_summ(n).ex_inbound_logic=[];
            
            Ripdis_summ(n).expre_inbound_perfvsday=[];
            Ripdis_summ(n).expre_inbound_rawperfvsday=rawperfvsday;
        end
    end
    
    
    % 3) Description
    Ripdis_summ(n).prefix=currprefix;
    
    if ~isempty(intersect(n,expidx))
        Ripdis_summ(n).group = 'Exp';
    end
    
    if ~isempty(intersect(n,conidx))
        Ripdis_summ(n).group = 'Con';
    end
    
    if ~isempty(intersect(n,normidx))
        Ripdis_summ(n).group = 'Nor';
    end
    
    % 4) Inbound Pers Error
    %Rip_inbound_summary(n).proportion_correct = frac_corr;
    %Rip_inbound_summary(n).proportion_perserr = frac_perserr;
    %Rip_inbound_summary(n).proportion_turnerr = frac_turnerr;
    %Rip_inbound_summary(n).number_trials = ntrajs_perday;
    
    
end


%******************************************************************
%% Separate by group and average
expcnt=0; concnt=0; norcnt=0;
expcnt_ex=0; concnt_ex=0;

for n=1:length(Ripdis_summ)
    
    currprefix=Ripdis_summ(n).prefix
    
    % Pad some points to get to 200
    if strcmp(Ripdis_summ(n).prefix,'REd')
        Ripdis_summ(n).out_smoothperf = [Ripdis_summ(n).out_smoothperf; zeros(13,1)];
        Ripdis_summ(n).out_curve = [Ripdis_summ(n).out_curve; zeros(13,1)];
        
        Ripdis_summ(n).out_smoothperf(188:200)=Ripdis_summ(n).out_smoothperf(175:187);
        Ripdis_summ(n).out_curve(198:200)=Ripdis_summ(n).out_curve(190:192);
        
    end
    
    if strcmp(Ripdis_summ(n).prefix,'REg')
       Ripdis_summ(n).out_smoothperf = [Ripdis_summ(n).out_smoothperf; zeros(25,1)];
       Ripdis_summ(n).out_curve = [Ripdis_summ(n).out_curve; zeros(20,1)];
       
       Ripdis_summ(n).out_smoothperf(176:200)=Ripdis_summ(n).out_smoothperf(151:175);
       Ripdis_summ(n).out_curve(181:200)=Ripdis_summ(n).out_curve(161:180);    
    end
    
    if strcmp(Ripdis_summ(n).prefix,'RNb')
        % keyboard;
        Ripdis_summ(n).out_smoothperf = [Ripdis_summ(n).out_smoothperf; zeros(18,1)];
        Ripdis_summ(n).out_curve = [Ripdis_summ(n).out_curve; zeros(8,1)];
        
        Ripdis_summ(n).out_smoothperf(183:200)=Ripdis_summ(n).out_smoothperf(161:178);
        Ripdis_summ(n).out_curve(192:200)=Ripdis_summ(n).out_curve(182:190);
        
    end
    
    switch Ripdis_summ(n).group
        
        
        
        case 'Exp'
            expcnt=expcnt+1;
            
            
            %             if n==2,
            %                 keyboard;
            %             end
            out_smoothperf_expall(expcnt,:)=Ripdis_summ(n).out_smoothperf(1:Ntraj_cutoff);
            out_curve_expall(expcnt,:)=Ripdis_summ(n).out_curve(1:Ntraj_cutoff);
            out_exp_lt(expcnt)=Ripdis_summ(n).out_learning_trial;
            out_exp_ls(expcnt)=Ripdis_summ(n).out_learning_sess;
            out_exp_ld(expcnt)=Ripdis_summ(n).out_learning_day;
            out_ntrajspd_expall(expcnt,:)=Ripdis_summ(n).out_ntrajs_perday(1:day_cutoff);
            out_ntrajvsd_expall(expcnt,:)=Ripdis_summ(n).out_ntraj_axis_day(1:day_cutoff);
            out_perfday_expall{expcnt}=Ripdis_summ(n).outbound_perf_day;
            out_rawperfvsday_expall{expcnt}=Ripdis_summ(n).out_rawperfvsday;
            
            in_smoothperf_expall(expcnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_expall(expcnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_exp_lt(expcnt)=Ripdis_summ(n).in_learning_trial;
            in_exp_ls(expcnt)=Ripdis_summ(n).in_learning_sess;
            in_exp_ld(expcnt)=Ripdis_summ(n).in_learning_day;
            in_ntrajspd_expall(expcnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_expall(expcnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_perfday_expall{expcnt}=Ripdis_summ(n).inbound_perf_day;
            in_rawperfvsday_expall{expcnt}=Ripdis_summ(n).in_rawperfvsday;
            
            % Last 2 days: Ncorrect and Ntrials
            out_last2d_logic_expall{expcnt}=Ripdis_summ(n).outbound_logic( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: Ripdis_summ(n).out_ntraj_axis_day(end));
            in_last2d_logic_expall{expcnt}=Ripdis_summ(n).inbound_logic( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: Ripdis_summ(n).in_ntraj_axis_day(end));
            
            if ~isempty(Ripdis_summ(n).ex_outbound_perf_day)
                expcnt_ex=expcnt_ex+1;
                
                out_ex_perfday_expall(expcnt_ex,:)=Ripdis_summ(n).ex_outbound_perf_day;
                out_ex_curve_expall{expcnt_ex}=Ripdis_summ(n).ex_out_curve;
                out_ex_logic_expall{expcnt_ex}=Ripdis_summ(n).ex_outbound_logic;
                out_expre_logic_expall{expcnt_ex}=Ripdis_summ(n).outbound_logic( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: Ripdis_summ(n).out_ntraj_axis_day(end));
                out_expre_rawperfvsday_expall{expcnt_ex} = Ripdis_summ(n).expre_outbound_rawperfvsday;
                out_expre_perfvsday_expall{expcnt_ex} = Ripdis_summ(n).expre_outbound_perfvsday;
                
                in_ex_perfday_expall(expcnt_ex,:)=Ripdis_summ(n).ex_inbound_perf_day;
                in_ex_curve_expall{expcnt_ex}=Ripdis_summ(n).ex_in_curve;
                in_ex_logic_expall{expcnt_ex}=Ripdis_summ(n).ex_inbound_logic;
                in_expre_logic_expall{expcnt_ex}=Ripdis_summ(n).inbound_logic( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: Ripdis_summ(n).in_ntraj_axis_day(end));
                in_expre_rawperfvsday_expall{expcnt_ex}  = Ripdis_summ(n).expre_inbound_rawperfvsday;
                in_expre_perfvsday_expall{expcnt_ex} = Ripdis_summ(n).expre_inbound_perfvsday;
                
            end
            
        case 'Con'
            concnt=concnt+1;
            
            out_smoothperf_conall(concnt,:)=Ripdis_summ(n).out_smoothperf(1:Ntraj_cutoff);
            out_curve_conall(concnt,:)=Ripdis_summ(n).out_curve(1:Ntraj_cutoff);
            out_con_lt(concnt)=Ripdis_summ(n).out_learning_trial;
            out_con_ls(concnt)=Ripdis_summ(n).out_learning_sess;
            out_con_ld(concnt)=Ripdis_summ(n).out_learning_day;
            out_ntrajspd_conall(concnt,:)=Ripdis_summ(n).out_ntrajs_perday(1:day_cutoff);
            out_ntrajvsd_conall(concnt,:)=Ripdis_summ(n).out_ntraj_axis_day(1:day_cutoff);
            out_perfday_conall{concnt}=Ripdis_summ(n).outbound_perf_day;
            out_rawperfvsday_conall{concnt}=Ripdis_summ(n).out_rawperfvsday;
            
            in_smoothperf_conall(concnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_conall(concnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_con_lt(concnt)=Ripdis_summ(n).in_learning_trial;
            in_con_ls(concnt)=Ripdis_summ(n).in_learning_sess;
            in_con_ld(concnt)=Ripdis_summ(n).in_learning_day;
            in_ntrajspd_conall(concnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_conall(concnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_perfday_conall{concnt}=Ripdis_summ(n).inbound_perf_day;
            in_rawperfvsday_conall{concnt}=Ripdis_summ(n).in_rawperfvsday;
            
             % Last 2 days: Ncorrect and Ntrials
            out_last2d_logic_conall{concnt}=Ripdis_summ(n).outbound_logic( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: Ripdis_summ(n).out_ntraj_axis_day(end));
            in_last2d_logic_conall{concnt}=Ripdis_summ(n).inbound_logic( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: Ripdis_summ(n).in_ntraj_axis_day(end));
            
            if ~isempty(Ripdis_summ(n).ex_outbound_perf_day)
                concnt_ex=concnt_ex+1;
                
                out_ex_perfday_conall(concnt_ex,:)=Ripdis_summ(n).ex_outbound_perf_day;
                out_ex_curve_conall{concnt_ex}=Ripdis_summ(n).ex_out_curve;
                out_ex_logic_conall{concnt_ex}=Ripdis_summ(n).ex_outbound_logic;
                out_expre_logic_conall{concnt_ex}=Ripdis_summ(n).outbound_logic( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: Ripdis_summ(n).out_ntraj_axis_day(end));
                out_expre_rawperfvsday_conall{concnt_ex}  = Ripdis_summ(n).expre_outbound_rawperfvsday;
                out_expre_perfvsday_conall{concnt_ex} = Ripdis_summ(n).expre_outbound_perfvsday;
                
                in_ex_perfday_conall(concnt_ex,:)=Ripdis_summ(n).ex_inbound_perf_day;
                in_ex_curve_conall{concnt_ex}=Ripdis_summ(n).ex_in_curve;
                in_ex_logic_conall{concnt_ex}=Ripdis_summ(n).ex_inbound_logic;
                in_expre_logic_conall{concnt_ex}=Ripdis_summ(n).inbound_logic( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: Ripdis_summ(n).in_ntraj_axis_day(end));
                in_expre_rawperfvsday_conall{concnt_ex}  = Ripdis_summ(n).expre_inbound_rawperfvsday;
                in_expre_perfvsday_conall{concnt_ex} = Ripdis_summ(n).expre_inbound_perfvsday;
                
            end
            
        case 'Nor'
            norcnt=norcnt+1;
            
            out_smoothperf_norall(norcnt,:)=Ripdis_summ(n).out_smoothperf(1:Ntraj_cutoff);
            out_curve_norall(norcnt,:)=Ripdis_summ(n).out_curve(1:Ntraj_cutoff);
            out_nor_lt(norcnt)=Ripdis_summ(n).out_learning_trial;
            out_nor_ls(norcnt)=Ripdis_summ(n).out_learning_sess;
            out_nor_ld(norcnt)=Ripdis_summ(n).out_learning_day;
            out_ntrajspd_norall(norcnt,:)=Ripdis_summ(n).out_ntrajs_perday(1:day_cutoff);
            out_ntrajvsd_norall(norcnt,:)=Ripdis_summ(n).out_ntraj_axis_day(1:day_cutoff);
            out_perfday_norall{norcnt}=Ripdis_summ(n).outbound_perf_day;
            out_rawperfvsday_norall{norcnt}=Ripdis_summ(n).out_rawperfvsday;
            
            in_smoothperf_norall(norcnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_norall(norcnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_nor_lt(norcnt)=Ripdis_summ(n).in_learning_trial;
            in_nor_ls(norcnt)=Ripdis_summ(n).in_learning_sess;
            in_nor_ld(norcnt)=Ripdis_summ(n).in_learning_day;
            in_ntrajspd_norall(norcnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_norall(norcnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_perfday_norall{norcnt}=Ripdis_summ(n).inbound_perf_day;
            in_rawperfvsday_norall{norcnt}=Ripdis_summ(n).in_rawperfvsday;
            
             % Last 2 days: Ncorrect and Ntrials
            out_last2d_logic_norall{norcnt}=Ripdis_summ(n).outbound_logic( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: Ripdis_summ(n).out_ntraj_axis_day(end));
            in_last2d_logic_norall{norcnt}=Ripdis_summ(n).inbound_logic( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: Ripdis_summ(n).in_ntraj_axis_day(end));
    end
    
end

%% Get averages

%**************************************************
% A.) Outbound
%**************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experimental Group

% Smooth Perf - Fraction Correct
out_smoothperf_exp = mean( out_smoothperf_expall,1);
out_smoothperf_experr = std( out_smoothperf_expall,1);
out_smoothperf_explowerr = out_smoothperf_exp-out_smoothperf_experr;
out_smoothperf_expuperr = out_smoothperf_exp+out_smoothperf_experr;
% Last 100 trials for Smooth Perf
out_smoothperfend_exp = out_smoothperf_exp(end-99:end);
out_smoothperfend_explowerr = out_smoothperf_explowerr(end-99:end);
out_smoothperfend_expuperr = out_smoothperf_expuperr(end-99:end);

% Performance Curves
out_curve_exp = mean( out_curve_expall,1);
out_curve_experr = std( out_curve_expall,1);
out_curve_explowerr = out_curve_exp-out_curve_experr;
out_curve_expuperr = out_curve_exp+out_curve_experr;
% Last 100 trials
out_curveend_exp = out_curve_exp(end-99:end);
out_curveend_explowerr = out_curve_explowerr(end-99:end);
out_curveend_expuperr = out_curve_expuperr(end-99:end);


% Ntrajsperday
out_ntrajspd_exp = mean( out_ntrajspd_expall,1);
out_ntrajspd_experr = std( out_ntrajspd_expall,1);
out_ntrajspd_explowerr = out_ntrajspd_exp-out_ntrajspd_experr;
out_ntrajspd_expuperr = out_ntrajspd_exp+out_ntrajspd_experr;

% Ntrajvsday
out_ntrajvsd_exp = mean( out_ntrajvsd_expall,1);
out_ntrajvsd_experr = std( out_ntrajvsd_expall,1);
out_ntrajvsd_explowerr = out_ntrajvsd_exp-out_ntrajvsd_experr;
out_ntrajvsd_expuperr = out_ntrajvsd_exp+out_ntrajvsd_experr;

% Performance on last 2 days
for i=1:expcnt
    out_perfday_expall{expcnt};
    out_perfday_exp(i) = mean(out_perfday_expall{i}(end-lastdays+1:end));
    out_rawperfday_exp(i) = mean(out_rawperfvsday_expall{i}(end-lastdays+1:end));
    out_last2d_Ntrials_exp(i) = length(out_last2d_logic_expall{i});
    out_last2d_Ncorrtrials_exp(i) = length(find(out_last2d_logic_expall{i}==1));
end

%Raw Performance vs day
% If all animals have upto day_cutoff days, simply Use this:
for i=1:expcnt
    out_rawperfvsday_expmat(i,:)=out_rawperfvsday_expall{i}(1:day_cutoff);
end

% Performance on extra days and Pre-extra days Limited to animals with
% Extra Days
for i=1:expcnt_ex
    ex_out_perfday_exp(i) = mean(out_ex_perfday_expall(i,:));
    expre_out_perfday_exp(i) = mean(out_expre_perfvsday_expall{i}(end-lastdays+1:end)); % Same as out_rawperfday_exp
    ex_out_rawperfday_exp(i) = mean(out_ex_logic_expall{i});
    expre_out_rawperfday_exp(i) = mean(out_expre_rawperfvsday_expall{i}(end-lastdays+1:end));
    
    % Number of correct trials and no of total trials for 2 extra and 2 pre-extra days
    ex_out_Ntrials_exp(i) = length(out_ex_logic_expall{i});
    ex_out_Ncorrtrials_exp(i) = length(find(out_ex_logic_expall{i}==1));
    expre_out_Ntrials_exp(i) = length(out_expre_logic_expall{i});
    expre_out_Ncorrtrials_exp(i) = length(find(out_expre_logic_expall{i}==1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control Group
out_smoothperf_con = mean(out_smoothperf_conall,1);
out_smoothperf_conerr = std(out_smoothperf_conall,1);
out_smoothperf_conlowerr = out_smoothperf_con-out_smoothperf_conerr;
out_smoothperf_conuperr = out_smoothperf_con+out_smoothperf_conerr;

% Last 100 trials for Smooth Perf
out_smoothperfend_con = out_smoothperf_con(end-99:end);
out_smoothperfend_conlowerr = out_smoothperf_conlowerr(end-99:end);
out_smoothperfend_conuperr = out_smoothperf_conuperr(end-99:end);

% Performance Curves
out_curve_con = mean(out_curve_conall,1);
out_curve_conerr = std(out_curve_conall,1);
out_curve_conlowerr = out_curve_con-out_curve_conerr;
out_curve_conuperr = out_curve_con+out_curve_conerr;
% Last 100 trials
out_curveend_con = out_curve_con(end-99:end);
out_curveend_conlowerr = out_curve_conlowerr(end-99:end);
out_curveend_conuperr = out_curve_conuperr(end-99:end);

% Ntrajsperday
out_ntrajspd_con = mean( out_ntrajspd_conall,1);
out_ntrajspd_conerr = std( out_ntrajspd_conall,1);
out_ntrajspd_conlowerr = out_ntrajspd_con-out_ntrajspd_conerr;
out_ntrajspd_conuperr = out_ntrajspd_con+out_ntrajspd_conerr;

% Ntrajvsday
out_ntrajvsd_con = mean( out_ntrajvsd_conall,1);
out_ntrajvsd_conerr = std( out_ntrajvsd_conall,1);
out_ntrajvsd_conlowerr = out_ntrajvsd_con-out_ntrajvsd_conerr;
out_ntrajvsd_conuperr = out_ntrajvsd_con+out_ntrajvsd_conerr;

% Performance on last 2 days
for i=1:concnt
    out_perfday_con(i) = mean(out_perfday_conall{i}(end-lastdays+1:end));
    out_rawperfday_con(i) = mean(out_rawperfvsday_conall{i}(end-lastdays+1:end));
    out_last2d_Ntrials_con(i) = length(out_last2d_logic_conall{i});
    out_last2d_Ncorrtrials_con(i) = length(find(out_last2d_logic_conall{i}==1));
end

%Raw Performance vs day
% If all animals have upto day_cutoff days, simply use this:
for i=1:concnt
    out_rawperfvsday_conmat(i,:)=out_rawperfvsday_conall{i}(1:day_cutoff);
end

% Performance on extra days and Pre-extra days Limited to animals with
% Extra Days
for i=1:concnt_ex
    ex_out_perfday_con(i) = mean(out_ex_perfday_conall(i,:));
    expre_out_perfday_con(i) = mean(out_expre_perfvsday_conall{i}(end-lastdays+1:end));
    ex_out_rawperfday_con(i) = mean(out_ex_logic_conall{i});
    expre_out_rawperfday_con(i) = mean(out_expre_rawperfvsday_conall{i}(end-lastdays+1:end));
    
    % Number of correct trials and no of total trials for 2 extra and 2 pre-extra days
    ex_out_Ntrials_con(i) = length(out_ex_logic_conall{i});
    ex_out_Ncorrtrials_con(i) = length(find(out_ex_logic_conall{i}==1));
    expre_out_Ntrials_con(i) = length(out_expre_logic_conall{i});
    expre_out_Ncorrtrials_con(i) = length(find(out_expre_logic_conall{i}==1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normal Group

% Smooth Perf - Fraction Correct
out_smoothperf_nor = mean( out_smoothperf_norall,1);
out_smoothperf_norerr = std( out_smoothperf_norall,1);
out_smoothperf_norlowerr = out_smoothperf_nor-out_smoothperf_norerr;
out_smoothperf_noruperr = out_smoothperf_nor+out_smoothperf_norerr;
% Last 100 trials for Smooth Perf
out_smoothperfend_nor = out_smoothperf_nor(end-99:end);
out_smoothperfend_norlowerr = out_smoothperf_norlowerr(end-99:end);
out_smoothperfend_noruperr = out_smoothperf_noruperr(end-99:end);

% Performance Curves
out_curve_nor = mean( out_curve_norall,1);
out_curve_norerr = std( out_curve_norall,1);
out_curve_norlowerr = out_curve_nor-out_curve_norerr;
out_curve_noruperr = out_curve_nor+out_curve_norerr;
% Last 100 trials
out_curveend_nor = out_curve_nor(end-99:end);
out_curveend_norlowerr = out_curve_norlowerr(end-99:end);
out_curveend_noruperr = out_curve_noruperr(end-99:end);


%Raw Performance vs day
% If all animals have upto day_cutoff days, simply use this:
for i=1:norcnt
    out_rawperfvsday_normat(i,:)=out_rawperfvsday_norall{i}(1:day_cutoff);
end
% % Method 2
% for d=1:day_cutoff
%     cnt=0;
%     for i=1:norcnt
%         if length(out_rawperfvsday_norall{i})>=d
%             cnt=cnt+1;
%             out_rawperfvsaday_norstr{d}(cnt)=out_rawperfvsday_norall{i}(d);
%         end
%     end
% end

% Performance on last 2 days
for i=1:norcnt
    out_perfday_nor(i) = mean(out_perfday_norall{i}(end-lastdays+1:end));
    out_rawperfday_nor(i) = mean(out_rawperfvsday_norall{i}(end-lastdays+1:end));
    out_last2d_Ntrials_nor(i) = length(out_last2d_logic_norall{i});
    out_last2d_Ncorrtrials_nor(i) = length(find(out_last2d_logic_norall{i}==1));
end

% Ntrajsperday
out_ntrajspd_nor = mean( out_ntrajspd_norall,1);
out_ntrajspd_norerr = std( out_ntrajspd_norall,1);
out_ntrajspd_norlowerr = out_ntrajspd_nor-out_ntrajspd_norerr;
out_ntrajspd_noruperr = out_ntrajspd_nor+out_ntrajspd_norerr;

% Ntrajvsday
out_ntrajvsd_nor = mean( out_ntrajvsd_norall,1);
out_ntrajvsd_norerr = std( out_ntrajvsd_norall,1);
out_ntrajvsd_norlowerr = out_ntrajvsd_nor-out_ntrajvsd_norerr;
out_ntrajvsd_noruperr = out_ntrajvsd_nor+out_ntrajvsd_norerr;



%**************************************************
% B) Inbound
%**************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experimental Group

% Smooth Perf - Fraction Correct
in_smoothperf_exp = mean( in_smoothperf_expall,1);
in_smoothperf_experr = std( in_smoothperf_expall,1);
in_smoothperf_explowerr = in_smoothperf_exp-in_smoothperf_experr;
in_smoothperf_expuperr = in_smoothperf_exp+in_smoothperf_experr;
% Last 100 trials for Smooth Perf
in_smoothperfend_exp = in_smoothperf_exp(end-99:end);
in_smoothperfend_explowerr = in_smoothperf_explowerr(end-99:end);
in_smoothperfend_expuperr = in_smoothperf_expuperr(end-99:end);


% Performance Curves
in_curve_exp = mean( in_curve_expall,1);
in_curve_experr = std( in_curve_expall,1);
in_curve_explowerr = in_curve_exp-in_curve_experr;
in_curve_expuperr = in_curve_exp+in_curve_experr;
% Last 100 trials
in_curveend_exp = in_curve_exp(end-99:end);
in_curveend_explowerr = in_curve_explowerr(end-99:end);
in_curveend_expuperr = in_curve_expuperr(end-99:end);

% Ntrajsperday

in_ntrajspd_exp = mean( in_ntrajspd_expall,1);
in_ntrajspd_experr = std( in_ntrajspd_expall,1);
in_ntrajspd_explowerr = in_ntrajspd_exp-in_ntrajspd_experr;
in_ntrajspd_expuperr = in_ntrajspd_exp+in_ntrajspd_experr;

% Ntrajvsday
in_ntrajvsd_exp = mean( in_ntrajvsd_expall,1);
in_ntrajvsd_experr = std( in_ntrajvsd_expall,1);
in_ntrajvsd_explowerr = in_ntrajvsd_exp-in_ntrajvsd_experr;
in_ntrajvsd_expuperr = in_ntrajvsd_exp+in_ntrajvsd_experr;

% Performance on last 2 days
for i=1:expcnt
    in_perfday_exp(i) = mean(in_perfday_expall{i}(end-lastdays+1:end));
    in_rawperfday_exp(i) = mean(in_rawperfvsday_expall{i}(end-lastdays+1:end));
    in_last2d_Ntrials_exp(i) = length(in_last2d_logic_expall{i});
    in_last2d_Ncorrtrials_exp(i) = length(find(in_last2d_logic_expall{i}==1));
end

%Raw Performance vs day
% If all animals have upto day_cutoff days, simply use this:
for i=1:expcnt
    in_rawperfvsday_expmat(i,:)=in_rawperfvsday_expall{i}(1:day_cutoff);
end

% Performance on extra days and Pre-extra days Limited to animals with
% Extra Days
for i=1:expcnt_ex
    ex_in_perfday_exp(i) = mean(in_ex_perfday_expall(i,:));
    expre_in_perfday_exp(i) = mean(in_expre_perfvsday_expall{i}(end-lastdays+1:end));
    ex_in_rawperfday_exp(i) = mean(in_ex_logic_expall{i});
    expre_in_rawperfday_exp(i) = mean(in_expre_rawperfvsday_expall{i}(end-lastdays+1:end));
    
    % Number of correct trials and no of total trials for 2 extra and 2 pre-extra days
    ex_in_Ntrials_con(i) = length(in_ex_logic_conall{i});
    ex_in_Ncorrtrials_con(i) = length(find(in_ex_logic_conall{i}==1));
    expre_in_Ntrials_con(i) = length(in_expre_logic_conall{i});
    expre_in_Ncorrtrials_con(i) = length(find(in_expre_logic_conall{i}==1));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control Group

in_smoothperf_con = mean( in_smoothperf_conall,1);
in_smoothperf_conerr = std( in_smoothperf_conall,1);
in_smoothperf_conlowerr = in_smoothperf_con-in_smoothperf_conerr;
in_smoothperf_conuperr = in_smoothperf_con+in_smoothperf_conerr;
% Last 100 trials for Smooth Perf
in_smoothperfend_con = in_smoothperf_con(end-99:end);
in_smoothperfend_conlowerr = in_smoothperf_conlowerr(end-99:end);
in_smoothperfend_conuperr = in_smoothperf_conuperr(end-99:end);

in_curve_con = mean( in_curve_conall,1);
in_curve_conerr = std( in_curve_conall,1);
in_curve_conlowerr = in_curve_con-in_curve_conerr;
in_curve_conuperr = in_curve_con+in_curve_conerr;
% Last 100 trials
in_curveend_con = in_curve_con(end-99:end);
in_curveend_conlowerr = in_curve_conlowerr(end-99:end);
in_curveend_conuperr = in_curve_conuperr(end-99:end);

in_ntrajspd_con = mean( in_ntrajspd_conall,1);
in_ntrajspd_conerr = std( in_ntrajspd_conall,1);
in_ntrajspd_conlowerr = in_ntrajspd_con-in_ntrajspd_conerr;
in_ntrajspd_conuperr = in_ntrajspd_con+in_ntrajspd_conerr;

in_ntrajvsd_con = mean( in_ntrajvsd_conall,1);
in_ntrajvsd_conerr = std( in_ntrajvsd_conall,1);
in_ntrajvsd_conlowerr = in_ntrajvsd_con-in_ntrajvsd_conerr;
in_ntrajvsd_conuperr = in_ntrajvsd_con+in_ntrajvsd_conerr;

% Performance on last 2 days
for i=1:concnt
    in_perfday_con(i) = mean(in_perfday_conall{i}(end-lastdays+1:end));
    in_rawperfday_con(i) = mean(in_rawperfvsday_conall{i}(end-lastdays+1:end));
    in_last2d_Ntrials_con(i) = length(in_last2d_logic_conall{i});
    in_last2d_Ncorrtrials_con(i) = length(find(in_last2d_logic_conall{i}==1));
end

%Raw Performance vs day
% If all animals have upto day_cutoff days, simply use this:
for i=1:concnt
    in_rawperfvsday_conmat(i,:)=in_rawperfvsday_conall{i}(1:day_cutoff);
end

% Performance on extra days and Pre-extra days Limited to animals with
% Extra Days
for i=1:concnt_ex
    ex_in_perfday_con(i) = mean(in_ex_perfday_conall(i,:));
    expre_in_perfday_con(i) = mean(in_expre_perfvsday_conall{i}(end-lastdays+1:end));
    ex_in_rawperfday_con(i) = mean(in_ex_logic_conall{i});
    expre_in_rawperfday_con(i) = mean(in_expre_rawperfvsday_conall{i}(end-lastdays+1:end));
    
    % Number of correct trials and no of total trials for 2 extra and 2 pre-extra days
    ex_in_Ntrials_con(i) = length(in_ex_logic_conall{i});
    ex_in_Ncorrtrials_con(i) = length(find(in_ex_logic_conall{i}==1));
    expre_in_Ntrials_con(i) = length(in_expre_logic_conall{i});
    expre_in_Ncorrtrials_con(i) = length(find(in_expre_logic_conall{i}==1));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normal Group

% Smooth Perf - Fraction Correct
in_smoothperf_nor = mean( in_smoothperf_norall,1);
in_smoothperf_norerr = std( in_smoothperf_norall,1);
in_smoothperf_norlowerr = in_smoothperf_nor-in_smoothperf_norerr;
in_smoothperf_noruperr = in_smoothperf_nor+in_smoothperf_norerr;
% Last 100 trials for Smooth Perf
in_smoothperfend_nor = in_smoothperf_nor(end-99:end);
in_smoothperfend_norlowerr = in_smoothperf_norlowerr(end-99:end);
in_smoothperfend_noruperr = in_smoothperf_noruperr(end-99:end);

% Performance Curves
in_curve_nor = mean( in_curve_norall,1);
in_curve_norerr = std( in_curve_norall,1);
in_curve_norlowerr = in_curve_nor-in_curve_norerr;
in_curve_noruperr = in_curve_nor+in_curve_norerr;
% Last 100 trials
in_curveend_nor = in_curve_nor(end-99:end);
in_curveend_norlowerr = in_curve_norlowerr(end-99:end);
in_curveend_noruperr = in_curve_noruperr(end-99:end);


% If all animals have upto day_cutoff days, simply use this:
for i=1:norcnt
    in_rawperfvsday_normat(i,:)=in_rawperfvsday_norall{i}(1:day_cutoff);
end
% % Method2
% for d=1:day_cutoff
%     cnt=0;
%     for i=1:norcnt
%         if length(in_rawperfvsday_norall{i})>=d
%             cnt=cnt+1;
%             in_rawperfvsaday_norstr{d}(cnt)=in_rawperfvsday_norall{i}(d);
%         end
%     end
% end

% Performance on last 2 days
for i=1:norcnt
    in_perfday_nor(i) = mean(in_perfday_norall{i}(end-lastdays+1:end));
    in_rawperfday_nor(i) = mean(in_rawperfvsday_norall{i}(end-lastdays+1:end));
    in_last2d_Ntrials_nor(i) = length(in_last2d_logic_norall{i});
    in_last2d_Ncorrtrials_nor(i) = length(find(in_last2d_logic_norall{i}==1));
end

% Ntrajsperday
in_ntrajspd_nor = mean( in_ntrajspd_norall,1);
in_ntrajspd_norerr = std( in_ntrajspd_norall,1);
in_ntrajspd_norlowerr = in_ntrajspd_nor-in_ntrajspd_norerr;
in_ntrajspd_noruperr = in_ntrajspd_nor+in_ntrajspd_norerr;

% Ntrajvsday
in_ntrajvsd_nor = mean( in_ntrajvsd_norall,1);
in_ntrajvsd_norerr = std( in_ntrajvsd_norall,1);
in_ntrajvsd_norlowerr = in_ntrajvsd_nor-in_ntrajvsd_norerr;
in_ntrajvsd_noruperr = in_ntrajvsd_nor+in_ntrajvsd_norerr;



%*******************************************************************
%*******************************************************************
%                          FIGURES
%*******************************************************************
%*******************************************************************

% Go to Figure Save Directory
cd(summdir);

if figcurve==1  % Choice For Plotting Fitted Curves
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1.) Outbound: Performance Curves
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    plot(out_curve_exp,'r','Linewidth',3);
    jbfill([1:Ntraj_cutoff],out_curve_expuperr,out_curve_explowerr,'r','r',1,1);
    plot(out_curve_con,'b','Linewidth',3);
    jbfill([1:Ntraj_cutoff],out_curve_conuperr,out_curve_conlowerr,'b','b',1,1);
    if ~isempty(normidx)
        plot(out_curve_nor,'k','Linewidth',3);
        jbfill([1:Ntraj_cutoff],out_curve_noruperr,out_curve_norlowerr,'k','k',1,1);
    end
    
    % Make Plot presentable
    title('Outbound Performance vs Trials','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0 Ntraj_cutoff 0 1.01]);
    plot([0:1:Ntraj_cutoff], 0.5*ones(size([0:1:Ntraj_cutoff])),'k--','Linewidth',1);
    set(gca,'TickDir','out');

    % Save
    if savefig1==1,
        figfile = [figdir,'0Outbound_Curve_Summ'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 2.) Outbound: Performance Curves - Last 100 trials
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    plot(out_curveend_exp,'r','Linewidth',3);
    jbfill([1:length(out_curveend_exp)],out_curveend_expuperr,out_curveend_explowerr,'r','r',1,1);
    plot(out_curveend_con,'b','Linewidth',3);
    jbfill([1:length(out_curveend_con)],out_curveend_conuperr,out_curveend_conlowerr,'b','b',1,1);
    if ~isempty(normidx)
        plot(out_curveend_nor,'k','Linewidth',3);
        jbfill(1:length(out_curveend_nor),out_curveend_noruperr,out_curveend_norlowerr,'k','k',1,1);
    end
    
    % Make Plot presentable
    title('Outbound Performance - Last 100 trials','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0 length(out_curveend_exp) 0 1.01]);
    plot([0:1:length(out_curveend_exp)], 0.5*ones(size([0:1:length(out_curveend_exp)])),'k--','Linewidth',1);
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0Outbound_CurveEnd_Summ'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % 3.) Inbound: Performance Curves
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    plot(in_curve_exp,'r','Linewidth',3);
    jbfill([1:Ntraj_cutoff],in_curve_expuperr,in_curve_explowerr,'r','r',1,1);
    plot(in_curve_con,'b','Linewidth',3);
    jbfill([1:Ntraj_cutoff],in_curve_conuperr,in_curve_conlowerr,'b','b',1,1);
    if ~isempty(normidx)
        plot(in_curve_nor,'k','Linewidth',3);
        jbfill([1:Ntraj_cutoff],in_curve_noruperr,in_curve_norlowerr,'k','k',1,1);
    end
    
    % Make Plot presentable
    title('Inbound Performance vs Trials','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0 Ntraj_cutoff 0 1.05]);
    plot([0:1:Ntraj_cutoff], 0.5*ones(size([0:1:Ntraj_cutoff])),'k--','Linewidth',1);
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0Inbound_Curve_Summ'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 4.) Inbound: Performance Curves - Last 100 trials
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    plot(in_curveend_exp,'r','Linewidth',3);
    jbfill(1:length(in_curveend_exp),in_curveend_expuperr,in_curveend_explowerr,'r','r',1,1);
    plot(in_curveend_con,'b','Linewidth',3);
    jbfill(1:length(in_curveend_con),in_curveend_conuperr,in_curveend_conlowerr,'b','b',1,1);
    if ~isempty(normidx)
        plot(in_curveend_nor,'k','Linewidth',3);
        jbfill(1:length(in_curveend_nor),in_curveend_noruperr,in_curveend_norlowerr,'k','k',1,1);
    end
    
    % Make Plot presentable
    title('Inbound Performance - Last 100 trials','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0 length(in_curveend_exp) 0 1.05]);
    plot([0:1:length(in_curveend_exp)], 0.5*ones(size([0:1:length(in_curveend_exp)])),'k--','Linewidth',1);
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0Inbound_CurveEnd_Summ'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
end % end if figcurve



%*******************************************************************
%*******************************************************************


if figsmooth==1
    
%     %1.) Outbound: Smooth Performance - Fraction Correct
%     figure; hold on;
%     if forppr==1
%         redimscreen_figforppr1;
%     else
%         redimscreen_figforppt1;
%     end
%     plot(out_smoothperf_exp,'r','Linewidth',3);
%     jbfill([1:Ntraj_cutoff],out_smoothperf_expuperr,out_smoothperf_explowerr,'r','r',1,1);
%     plot(out_smoothperf_con,'b','Linewidth',3);
%     jbfill([1:Ntraj_cutoff],out_smoothperf_conuperr,out_smoothperf_conlowerr,'b','b',1,1);
%     if ~isempty(normidx)
%         plot(out_smoothperf_nor,'k','Linewidth',3);
%         jbfill([1:Ntraj_cutoff],out_smoothperf_noruperr,out_smoothperf_norlowerr,'k','k',1,1);
%     end
%     
%     % Make Plot presentable
%     title('Outbound Smooth Performance vs Trials','FontSize',tfont,'Fontweight','normal');
%     xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
%     ylabel('Proportion correct','FontSize',yfont,'Fontweight','normal');
%     %legend('Location','NorthEast');
%     axis([0 Ntraj_cutoff 0 1.1]);
%     plot([0:1:Ntraj_cutoff], 0.5*ones(size([0:1:Ntraj_cutoff])),'k--','Linewidth',1);
%     set(gca,'TickDir','out');
%     
%     % Save
%     if savefig1==1,
%         figfile = [figdir,'Outbound_SmoothPerf_Summ'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     % 2.) Outbound: Smooth Performance Curves - Last 100 trials
%   
%     figure; hold on; 
%     if forppr==1
%         redimscreen_figforppr1;
%     else
%         redimscreen_figforppt1;
%     end
%     plot(out_smoothperfend_exp,'r','Linewidth',3);
%     jbfill([1:length(out_smoothperfend_exp)],out_smoothperfend_expuperr,out_smoothperfend_explowerr,'r','r',1,1);
%     plot(out_smoothperfend_con,'b','Linewidth',3);
%     jbfill([1:length(out_smoothperfend_con)],out_smoothperfend_conuperr,out_smoothperfend_conlowerr,'b','b',1,1);
%     
%     % Make Plot presentable
%     title('Outbound Smooth Performance - Last 100 trials','FontSize',tfont,'Fontweight','normal');
%     xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
%     ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
%     %legend('Location','NorthEast');
%     axis([0 length(out_smoothperfend_exp) 0 1.05]);
%     plot([0:1:length(out_smoothperfend_exp)], 0.5*ones(size([0:1:length(out_smoothperfend_exp)])),'k--','Linewidth',1);
%     set(gca,'TickDir','out');
%     
%     % Save
%     if savefig1==1,
%         figfile = [figdir,'Outbound_SmoothPerfEnd_Summ'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%     %3.) Inbound Smooth Performance
%     
%     figure; hold on; 
%     if forppr==1
%         redimscreen_figforppr1;
%     else
%         redimscreen_figforppt1;
%     end
%     plot(in_smoothperf_exp,'r','Linewidth',3);
%     jbfill([1:Ntraj_cutoff],in_smoothperf_expuperr,in_smoothperf_explowerr,'r','r',1,1);
%     plot(in_smoothperf_con,'b','Linewidth',3);
%     jbfill([1:Ntraj_cutoff],in_smoothperf_conuperr,in_smoothperf_conlowerr,'b','b',1,1);
%     
%     % Make Plot presentable
%     title('Inbound Smooth Performance vs Trials','FontSize',tfont,'Fontweight','normal');
%     xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
%     ylabel('Proportion correct','FontSize',yfont,'Fontweight','normal');
%     %legend('Location','NorthEast');
%     axis([0 Ntraj_cutoff 0 1.1]);
%     plot([0:1:Ntraj_cutoff], 0.5*ones(size([0:1:Ntraj_cutoff])),'k--','Linewidth',1);
%     set(gca,'TickDir','out');
%     
%     % Save
%     if savefig1==1,
%         figfile = [figdir,'Inbound_SmoothPerf_Summ'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     % 4.) Inbound: Smooth Performance Curves - Last 100 trials  
%     
%     figure; hold on; 
%     if forppr==1
%         redimscreen_figforppr1;
%     else
%         redimscreen_figforppt1;
%     end
%     plot(in_smoothperfend_exp,'r','Linewidth',3);
%     jbfill([1:length(in_smoothperfend_exp)],in_smoothperfend_expuperr,in_smoothperfend_explowerr,'r','r',1,1);
%     plot(in_smoothperfend_con,'b','Linewidth',3);
%     jbfill([1:length(in_smoothperfend_con)],in_smoothperfend_conuperr,in_smoothperfend_conlowerr,'b','b',1,1);
%     
%     % Make Plot presentable
%     title('Inbound Smooth Performance - Last 100 trials','FontSize',tfont,'Fontweight','normal');
%     xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
%     ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
%     %legend('Location','NorthEast');
%     axis([0 length(in_smoothperfend_exp) 0 1.05]);
%     plot([0:1:length(in_smoothperfend_exp)], 0.5*ones(size([0:1:length(in_smoothperfend_exp)])),'k--','Linewidth',1);
%     set(gca,'TickDir','out');
%     
%     % Save
%     if savefig1==1,
%         figfile = [figdir,'Inbound_SmoothPerf_Summ'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 5.) Outbound Raw Performance vs day for Exp, Con and Normal Groups
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    meanexp = mean(out_rawperfvsday_expmat,1);
    errexp = sem(out_rawperfvsday_expmat,1);
    meancon = mean(out_rawperfvsday_conmat,1);
    errcon = sem(out_rawperfvsday_conmat,1);
    
    errorbar([1.1:8.1],meanexp,errexp,'r.-','Linewidth',2,'MarkerSize',20);
    errorbar([1:8],meancon,errcon,'b.-','Linewidth',2,'MarkerSize',20);
    
    if ~isempty(normidx)
        meannor = mean(out_rawperfvsday_normat,1);
        errnor = sem(out_rawperfvsday_normat,1);
        errorbar([0.9:7.9],meannor,errnor,'k.-','Linewidth',2,'MarkerSize',20);
        %         for d=1:length(out_rawperfvsaday_norstr),
        %             normean(d) = mean(out_rawperfvsaday_norstr{d});
        %             norstd(d) = std(out_rawperfvsaday_norstr{d});
        %         end
        %         errorbar([0.9:7.9],normean,norstd,'k.-','Linewidth',2,'MarkerSize',18);
    end
    
    % Make Plot presentable
    title('Outbound Performance vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Proportion Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 8.5 0 1.0]);
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0Outbound_RawPerfVsDay_semerr'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %% Stats
    
     out_perfday_exp = out_rawperfvsday_expmat(:);
     out_perfday_con = out_rawperfvsday_conmat(:);
     out_perfday_nor = out_rawperfvsday_normat(:);
    
    % When I had 4 animals in each group
    % ----------------------------------
    
%     % Anova between mean curves
%     % Anova1 - Effect of group only.  
%      pout_expcon_rawperfvsday = anova1([meanexp',meancon'],[],'off'); % p = 0.0114
%      pout_expnor_rawperfvsday = anova1([meanexp',meannor'],[],'off'); % p = 0.0141;
%      pout_connor_rawperfvsday = anova1([meancon',meannor'],[],'off'); % p = 0.6345;
%      pout_expconnor_rawperfvsday = anova1([meanexp',meancon',meannor'],[],'off'); % p = 0.0128;
%     
%     % Anova with group X day effects. 4 reps=animals each for each group and day.
%     % First arrange matrices - Use out_rawperfvsday_expmat(4X8) or out_rawperfvsday_expall(1X4 cell)
%     % No need of looping - get it directly
%     %nan=4;
%     %for day = 1:size(out_rawperfvsday_expmat,2)
%     %    for an = 1:size(out_rawperfvsday_expmat,1)
%     %        expperfday(((day-1)*nan)+an) = out_rawperfvsday_expmat(an,day);
%     %    end
%     %end
%     
     
%     % anova 2
%     pout_expcon_grpday = anova2([out_perfday_exp,out_perfday_con],4,'off'); % pgrp=0.0000, pday=0.0002, pint=0.0146
%     pout_expnor_grpday = anova2([out_perfday_exp,out_perfday_nor],4,'off'); % pgrp=0.0000, pday=0.0000, pint=0.1558
%     pout_connor_grpday = anova2([out_perfday_con,out_perfday_nor],4,'off'); % pgrp=0.2935, pday=0.0000, pint=0.6276
%     [pout_expconnor_grpday, table, stats] = anova2([out_perfday_exp, out_perfday_con,out_perfday_nor],4,'off'); % pgrp=0.0000, pday=0.0000, pint=0.0499
%     
%     % multcompare on the anova2
%     % Using tukey-kramer = hsd by default
%     % 1a) - % compare grps: all days. Grp 2 and 3 diff from1 
%     [c,m,h,nms] = multcompare(stats,'estimate','column','display','off'); 
%     % 1b) - - % compare days: all grps. Days 3-8 diffeent from 1. Days 5,7,8 different from 1
%     [c,m,h,nms] = multcompare(stats,'estimate','row','display','off'); 
%     
%     % Using bonferroni = lsd. Same results 
%     % 1a) - % compare grps: all days. Grp 2 and 3 diff from1 
%     [c,m,h,nms] = multcompare(stats,'estimate','column','ctype','bonferroni','display','off'); 
%     % 1b) - - % compare days: all grps. Days 3-8 diffeent from 1. Days 5,7,8 different from 1
%     [c,m,h,nms] = multcompare(stats,'estimate','row','ctype','bonferroni','display','off'); 
    
    % ------------ End 4 animals in each group ----------------------------------------
    
    % anova repated measure
    pout_expcon_grpdayrm = anova_rm({out_rawperfvsday_expmat out_rawperfvsday_conmat},'off');
    pout_expnor_grpdayrm = anova_rm({out_rawperfvsday_expmat out_rawperfvsday_normat},'off');
    pout_connor_grpdayrm = anova_rm({out_rawperfvsday_normat out_rawperfvsday_conmat},'off'); %pgrp=0.4448; ptime=0; pint=0.5633
    pout_expconnor_grpdayrm = anova_rm({out_rawperfvsday_expmat out_rawperfvsday_conmat out_rawperfvsday_normat},'off'); %pgrp=0.0008, ptime=0; pint=0.0026;psub=0.0051  
    
    % Combine con and nor grps
    out_perfday = [out_perfday_exp;out_perfday_con;out_perfday_nor];
    grp = [0*ones(48,1);1*ones(64,1)]; % 6x8 Exp; 4x8 + 4x8 Con
    dayel_exp = [1*ones(6,1);2*ones(6,1);3*ones(6,1);4*ones(6,1);5*ones(6,1);6*ones(6,1);7*ones(6,1);8*ones(6,1)]; 
    dayel = [1*ones(4,1);2*ones(4,1);3*ones(4,1);4*ones(4,1);5*ones(4,1);6*ones(4,1);7*ones(4,1);8*ones(4,1)]; 
    dy = [dayel_exp; repmat(dayel,2,1)];
     % anova 2
    pout_expconandnor_grpday = anovan(out_perfday,{grp dy},'model','interaction','varnames',{'Grp';'Day'},'display','off'); % pgrp=0.0000, pday=0.0000, pint=0.0011
     % anova repated measure
    pout_expconandnor_grpdayrm = anova_rm({out_rawperfvsday_expmat [out_rawperfvsday_conmat;out_rawperfvsday_normat]},'off'); % pgrp=0.0001, pday=0.0000, pint=0.0002,psub=0.005
    
    % Pairwise t-test 
    for i=1:day_cutoff
        currexp = out_rawperfvsday_expmat(:,i);
        currcon = out_rawperfvsday_conmat(:,i);
        [p_rawperfvsday(i),h_rawperfvsday(i)]=ranksum(currexp,currcon);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 6.) Inbound Raw Performance vs day for Exp, Con and Normal Groups
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    meanexp = mean(in_rawperfvsday_expmat,1);
    errexp = sem(in_rawperfvsday_expmat,1);
    meancon = mean(in_rawperfvsday_conmat,1);
    errcon = sem(in_rawperfvsday_conmat,1);
    errorbar([1.1:8.1],meanexp,errexp,'r.-','Linewidth',2,'MarkerSize',20);
    errorbar([1:8],meancon,errcon,'b.-','Linewidth',2,'MarkerSize',20);
    
    if ~isempty(normidx)
        meannor = mean(in_rawperfvsday_normat,1);
        errnor = sem(in_rawperfvsday_normat,1);
        errorbar([0.9:7.9],meannor,errnor,'k.-','Linewidth',2,'MarkerSize',20);
        %         for d=1:length(in_rawperfvsaday_norstr),
        %             normean(d) = mean(in_rawperfvsaday_norstr{d});
        %             norstd(d) = std(in_rawperfvsaday_norstr{d});
        %         end
        %         errorbar([0.9:7.9],normean,norstd,'k.-','Linewidth',2,'MarkerSize',18);
    end
    
    % Make Plot presentable
    title('Inbound Performance vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Proportion Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 8.5 0 1.1]);
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0Inbound_RawPerfVsDay_semerr'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %% Stats
    
    in_perfday_exp = in_rawperfvsday_expmat(:);
    in_perfday_con = in_rawperfvsday_conmat(:);
    in_perfday_nor = in_rawperfvsday_normat(:);
    
    % When I had 4 animals in each group
    % ----------------------------------
%     
%     % Anova between mean curves
%     % Anova1 - Effect of group only.  
%     pin_expcon_rawperfvsday = anova1([meanexp',meancon'],[],'off'); % p = 0.5980
%     pin_expnor_rawperfvsday = anova1([meanexp',meannor'],[],'off'); % p = 0.9776;
%     pin_connor_rawperfvsday = anova1([meancon',meannor'],[],'off'); % p = 0.5395;
%     pin_expconnor_rawperfvsday = anova1([meanexp',meancon',meannor'],[],'off'); % p = 0.7973;
%     
%     % Anova with group X day effects. 4 reps=animals each for each group
%     % and day.
%     
%     % Anova2
%     pin_expcon_grpday = anova2([in_perfday_exp,in_perfday_con],4,'off'); % pgrp=0.0235, pday=0.0000, pint=0.9791
%     pin_expnor_grpday = anova2([in_perfday_exp,in_perfday_nor],4,'off'); % pgrp=0.9299, pday=0.0000, pint=0.7089
%     pin_connor_grpday = anova2([in_perfday_con,in_perfday_nor],4,'off'); % pgrp=0.0607, pday=0.0000, pint=0.8270
%     [pin_expconnor_grpday,table,stats] = anova2([in_perfday_exp, in_perfday_con,in_perfday_nor],4,'off'); % pgrp=0.0649, pday=0.0000, pint=0.9278
%     % multcompare on the anova2
%     % Using tukey-kramer = hsd by default
%     % 1a) - % compare grps: all days. All grps same
%     [c,m,h,nms] = multcompare(stats,'estimate','column','display','off'); 
%     % 1b) - - % compare days: all grps. Days 2-8 different from day1.  Days 4-8 different from day1. 
%     [c,m,h,nms] = multcompare(stats,'estimate','row','display','off'); 
%     
%     % Using bonferroni = lsd. Slightly different results for day2 comparisons
%     % 1a) - % compare grps: all days. All grps same
%     [c,m,h,nms] = multcompare(stats,'estimate','column','ctype','bonferroni','display','off'); 
%     % 1b) - - % compare days: all grps. Days 2-8 different from day1. Days 5,7,8 different from day1
%     [c,m,h,nms] = multcompare(stats,'estimate','row','ctype','bonferroni','display','off'); 
    
     % ------------ End 4 animals in each group ----------------------------------------
    
    % anova repated measure
    pin_expcon_grpdayrm = anova_rm({in_rawperfvsday_expmat in_rawperfvsday_conmat},'off'); % pgrp=0.2102, pday=0.0000, pint=0.8802
    pin_expnor_grpdayrm = anova_rm({in_rawperfvsday_expmat in_rawperfvsday_normat},'off');
    pin_connor_grpdayrm = anova_rm({in_rawperfvsday_normat in_rawperfvsday_conmat},'off'); % pgrp=0.1409, pday=0.0000, pint=0.7947
    pin_expconnor_grpdayrm = anova_rm({in_rawperfvsday_expmat in_rawperfvsday_conmat in_rawperfvsday_normat},'off'); % pgrp=0.3357, pday=0.0000, pint=0.5987
    
    % Combine con and nor grps
    in_perfday = [in_perfday_exp;in_perfday_con;in_perfday_nor];
    grp = [0*ones(48,1);1*ones(64,1)]; % 5x8 Exp; 4x8 + 4x8 Con
    dayel_exp = [1*ones(6,1);2*ones(6,1);3*ones(6,1);4*ones(6,1);5*ones(6,1);6*ones(6,1);7*ones(6,1);8*ones(6,1)]; 
    dayel = [1*ones(4,1);2*ones(4,1);3*ones(4,1);4*ones(4,1);5*ones(4,1);6*ones(4,1);7*ones(4,1);8*ones(4,1)]; 
    dy = [dayel_exp; repmat(dayel,2,1)];
    %anova2
    pin_expconandnor_grpday = anovan(in_perfday,{grp dy},'model','interaction','varnames',{'Grp';'Day'},'display','off'); % pgrp=0.199, pday=0.0000, pint=0.6678
    %anova_rm
    pin_expconandnor_grpdayrm = anova_rm({in_rawperfvsday_expmat [in_rawperfvsday_conmat;in_rawperfvsday_normat]},'off'); % pgrp=0.1817, pday=0.0000, pint=0.4788
    
    
end % end if figsmooth



%*******************************************************************
%*******************************************************************




if figout==1
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1.) Outbound: Learning Trial, Session and Day and Ntrials
    % Trial
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    lt(1,:)=[mean(out_exp_lt),mean(out_con_lt)];
    %bar(lt);
    plot(2,out_con_lt(1),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,out_con_lt(2),'bd','MarkerSize',16,'Linewidth',3);
    plot(1.9,out_con_lt(3),'bd','MarkerSize',16,'Linewidth',3);
    plot(2.1,out_con_lt(4),'bd','MarkerSize',16,'Linewidth',3);
    
    plot(1,out_exp_lt(1),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_exp_lt(2),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_exp_lt(3),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_exp_lt(4),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,out_exp_lt(5),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,out_exp_lt(6),'ro','MarkerSize',16,'Linewidth',3);
    
    %errorbar(0.8,mean(out_con_lt),sem(out_con_lt),'bd','MarkerSize',16,'Linewidth',4);
    %errorbar(1.2,mean(out_exp_lt),sem(out_exp_lt),'rd','MarkerSize',16,'Linewidth',4)
    
    axis([0 3 -5 200]);
    set(gca,'XTick',[1 2]);
    title('Outbound Learning Trial','FontSize',tfont,'Fontweight','normal');
    ylabel('Trials to reach Learning Criterion','FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');

    if ~isempty(normidx)
        for i=1:length(out_nor_lt)
            jitter(i)=0+0.05*(i-1);
        end
        %plot(3*ones(size(out_nor_lt))+jitter,out_nor_lt,'ks','MarkerSize',16,'Linewidth',4,'MarkerFaceColor','w');
        plot(2.9,out_nor_lt(1),'ks','MarkerSize',16,'Linewidth',3);
        plot(3,out_nor_lt(2),'ks','MarkerSize',16,'Linewidth',3);
        plot(3.1,out_nor_lt(3),'ks','MarkerSize',16,'Linewidth',3);
        plot(3,out_nor_lt(4),'ks','MarkerSize',16,'Linewidth',3);
        axis([0 4 -5 205]);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    end
    
    % Stats
    pout_expcon_lt = ranksum(out_con_lt,out_exp_lt); %p=0.0095
    pout_expnor_lt = ranksum(out_nor_lt,out_exp_lt); %p=0.0095
    pout_connor_lt = ranksum(out_con_lt,out_nor_lt); %p=0.09
    pout_expvsconnor_lt = ranksum(out_exp_lt,[out_con_lt,out_nor_lt]); %p=0.00066 (6.6e-4)
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0Outbound_LearningTrial'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 2.) Day/Sess
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    ls(1,:)=[mean(out_exp_ls),mean(out_con_ls)];
    ld(1,:)=[mean(out_exp_ld),mean(out_con_ld)];
    %bar(lt);
    plot(2,out_con_ld(1),'bd','MarkerSize',16,'Linewidth',3);
    plot(2.1,out_con_ld(2),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,out_con_ld(3),'bd','MarkerSize',16,'Linewidth',3);
    plot(1.9,out_con_ld(4),'bd','MarkerSize',16,'Linewidth',3);
    
    plot(1.1,out_exp_ld(1),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_exp_ld(2),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_exp_ld(3),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_exp_ld(4),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,out_exp_ld(5),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,out_exp_ld(6),'ro','MarkerSize',16,'Linewidth',3);
    %errorbar(0.8,mean(out_con_ld),sem(out_con_ld),'bd','MarkerSize',16,'Linewidth',2);
    %errorbar(1.2,mean(out_exp_ld),sem(out_exp_ld),'rd','MarkerSize',16,'Linewidth',2)
    
    axis([0 3 0.5 8.5]);
    set(gca,'XTick',[1 2]);
    title('Outbound Learning Day','FontSize',tfont,'Fontweight','normal');
    ylabel('Days to reach Learning Criterion','FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    if ~isempty(normidx)
        for i=1:length(out_nor_ld)
            jitter(i)=0+0.1*(i-1);
        end
        jitter=[-0.1 0 0.1 0.2];
        plot(3*ones(size(out_nor_ld))+jitter,out_nor_ld,'ks','MarkerSize',16,'Linewidth',3);
        axis([0 4 0.5 8.5]);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    end
    
    % Stats
    pout_expcon_ld = ranksum(out_con_ld,out_exp_ld); %p=0.0095
    pout_expnor_ld = ranksum(out_nor_ld,out_exp_ld); %p=0.0095
    pout_connor_ld = ranksum(out_con_ld,out_nor_ld); %p=1
    pout_expvsconnor_lt = ranksum(out_exp_ld,[out_con_ld,out_nor_ld]); %p=0.0006 (6.6e-4)

    % Save
    if savefig1==1,
        figfile = [figdir,'0Outbound_LearningDay'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 3.) Total Ntrials in day_cutoff days
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    out_TotalNtr_expall =sum(out_ntrajspd_expall,2);
    out_TotalNtr_conall =sum(out_ntrajspd_conall,2);
    out_TotalNtr_norall =sum(out_ntrajspd_norall,2);
    %plot(out_TotalNtr_conall,2*ones(size(out_TotalNtr_conall)),'bd','MarkerSize',16,'Linewidth',4,'MarkerFaceColor','b');
    plot(2,out_TotalNtr_conall(1),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,out_TotalNtr_conall(2),'bd','MarkerSize',16,'Linewidth',3);
    plot(2.1,out_TotalNtr_conall(3),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,out_TotalNtr_conall(4),'bd','MarkerSize',16,'Linewidth',3);
    %plot(out_TotalNtr_expall,1*ones(size(out_TotalNtr_expall)),'ro','MarkerSize',16,'Linewidth',4);
    plot(1.07,out_TotalNtr_expall(1),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_TotalNtr_expall(2),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_TotalNtr_expall(3),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,out_TotalNtr_expall(4),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,out_TotalNtr_expall(5),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,out_TotalNtr_expall(6),'ro','MarkerSize',16,'Linewidth',3);
    
    axis([0 3 0 550]);
    set(gca,'XTick',[1 2]);
    title(['Total Number of Outbound Trials'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Number of Trials'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    if ~isempty(normidx)
        for i=1:length(out_TotalNtr_norall)
            jitter(i)=0+0.05*(i-1);
        end
        plot(3*ones(size(out_TotalNtr_norall))+jitter',out_TotalNtr_norall,'ks','MarkerSize',16,'Linewidth',3);
        axis([0 4 0 550]);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    end
    
    % Stats
    pout_expcon_Ntrials = ranksum(out_TotalNtr_conall,out_TotalNtr_expall);%p=0.2571
    pout_expnor_Ntrials = ranksum(out_TotalNtr_norall,out_TotalNtr_expall);%p=0.3524
    pout_connor_Ntrials = ranksum(out_TotalNtr_conall,out_TotalNtr_norall);%p=1
    pout_expvsconnor_lt = ranksum(out_TotalNtr_expall,[out_TotalNtr_conall;out_TotalNtr_norall]);%p=0.1812
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0Outbound_NTrials'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 4.) Outbound Ntrials vs days
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    %     plot(out_ntrajspd_exp,'r.-','Linewidth',2,'MarkerSize',18);
    %     jbfill([1:length(out_ntrajvsd_exp)],out_ntrajvsd_expuperr,out_ntrajvsd_explowerr,'r','r',1,0.3);
    %     plot(out_ntrajspd_con,'b.-','Linewidth',2,'MarkerSize',18);
    %     jbfill([1:length(out_ntrajvsd_con)],out_ntrajvsd_conuperr,out_ntrajvsd_conlowerr,'b','b',1,0.3);
    errorbar([1:8],out_ntrajspd_exp,out_ntrajspd_experr,'ro-','Linewidth',2,'MarkerSize',12);
    errorbar([1.1:8.1],out_ntrajspd_con,out_ntrajspd_conerr,'bd-','Linewidth',2,'MarkerSize',12);
    [panovarm_Ntrdays_ConExp] = anova_rm({out_ntrajspd_expall out_ntrajspd_conall},'off');
    %ConExp: pgrp=0.2148, ptime=0, pint=0.4829
    [panovarm_Ntrdays_ConExpNor] = anova_rm({out_ntrajspd_expall out_ntrajspd_conall out_ntrajspd_norall},'off');
    %ConExpNor: pgrp=0.3948, ptime=0, pint=0.3119
    %if ~isempty(normidx)
    %    errorbar([0.9:7.9],out_ntrajspd_nor,out_ntrajspd_norerr,'k.-','Linewidth',2,'MarkerSize',20);
    %end
    
    % Make Plot presentable
    title('Number of Outbound Trials vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Days','FontSize',xfont,'Fontweight','normal');
    ylabel('Outbound Trials','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 8.5 0 80]);
    set(gca,'XTick',[1:8],'XTickLabel',{[1:8]'});
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0ConvsExp_Outbound_NTrialsPerDayVsDay'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 5.) Outbound Performance in last 2 days - USE RAW PERF ONLY
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    %     To Use Raw Performance rather than fitted performance, use out_rawperfday_con
    %     Raw Performance
        plot(2,out_rawperfday_con(1),'bd','MarkerSize',16,'Linewidth',3);
        plot(2.1,out_rawperfday_con(2),'bd','MarkerSize',16,'Linewidth',3);
        plot(2,out_rawperfday_con(3),'bd','MarkerSize',16,'Linewidth',3);
        plot(1.9,out_rawperfday_con(4),'bd','MarkerSize',16,'Linewidth',3);
        for i=1:length(out_rawperfday_exp)
            jitter(i)=0+0.05*(i-1);
        end
        %plot(1*ones(size(out_rawperfday_exp))+jitter,out_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',3);
        plot(1,out_rawperfday_exp(1),'ro','MarkerSize',16,'Linewidth',3);
        plot(1,out_rawperfday_exp(2),'ro','MarkerSize',16,'Linewidth',3);
        plot(0.85,out_rawperfday_exp(3),'ro','MarkerSize',16,'Linewidth',3);
        plot(1.0,out_rawperfday_exp(4),'ro','MarkerSize',16,'Linewidth',3);
        plot(1.2,out_rawperfday_exp(5),'ro','MarkerSize',16,'Linewidth',3);
        plot(1.2,out_rawperfday_exp(6),'ro','MarkerSize',16,'Linewidth',3);
    
      % Fitted Performance
%     plot(2,out_perfday_con(1),'bd','MarkerSize',16,'Linewidth',3);
%     plot(2.1,out_perfday_con(2),'bd','MarkerSize',16,'Linewidth',3);
%     plot(2,out_perfday_con(3),'bd','MarkerSize',16,'Linewidth',3);
%     plot(1.9,out_perfday_con(4),'bd','MarkerSize',16,'Linewidth',3);
%     plot(1*ones(size(out_perfday_expall)),out_perfday_exp,'ro','MarkerSize',16,'Linewidth',3);
    
    axis([0 3 0 1.05]);
    set(gca,'XTick',[1 2]);
    title(['Outbound Performance in last ' num2str(lastdays) ' days'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    if ~isempty(normidx)
        jitter=[];
        for i=1:length(out_rawperfday_nor)
            jitter(i)=0+0.05*(i-1);
        end
        % Raw Performance
        plot(3*ones(size(out_rawperfday_nor))+jitter,out_perfday_nor,'ks','MarkerSize',16,'Linewidth',3);
        % Fitted Performance
%       plot(3*ones(size(out_perfday_nor))+jitter,out_perfday_nor,'ks','MarkerSize',16,'Linewidth',3);
        
        axis([0 4 0 1.05]);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    end
    
    %Stats
    pout_expcon_last2days = ranksum(out_rawperfday_con,out_rawperfday_exp); %p=0.016
    pout_expnor_last2days = ranksum(out_rawperfday_nor,out_rawperfday_exp); %p=0.016
    pout_connor_last2days = ranksum(out_rawperfday_nor,out_rawperfday_con); %p=0.49
    pout_expvsconnor_last2days = ranksum(out_rawperfday_exp,[out_rawperfday_con, out_rawperfday_nor]); %p=0.0016

    % Save
    if savefig1==1,
        figfile = [figdir,'0OutboundRawPerf_Last2Days'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % 5b.) Outbound Performance in last 2 days - Plot Sem instead of points
    % for all 3 groups - use Raw Perf Only
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % To Use Raw Performance, use out_rawperfday_con
    errorbar(2,mean(out_rawperfday_con), sem(out_rawperfday_con),'bd','MarkerSize',16,'Linewidth',3);
    errorbar(1,mean(out_rawperfday_exp), sem(out_rawperfday_exp),'ro','MarkerSize',16,'Linewidth',3);
    errorbar(3,mean(out_rawperfday_nor), sem(out_rawperfday_nor),'ks','MarkerSize',16,'Linewidth',3);
 
    axis([0 4 0 1.05]);
    set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    title(['Outbound Performance in last ' num2str(lastdays) ' days'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');

    % Save
    if savefig1==1,
        figfile = [figdir,'OutboundRawPerf_Last2Days_witherr'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 5c.)  Prop Correct ztest Ratio: Outbound Performance in last 2 days - 
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
     
     % ([395 481] [261 371]) - old, for 4 anim each
     % 
    [pout_expcon_ratio, out_stderrp_expcon, out_stderrp_exp, out_stderrp_con]= sj_ztestprop2([sum(out_last2d_Ncorrtrials_con) sum(out_last2d_Ntrials_con)] , [sum(out_last2d_Ncorrtrials_exp) sum(out_last2d_Ntrials_exp)]);
    pout_expcon_ratio; % 0.0014
    out_stderrp_expcon; % 0.0358
    
     % ([321 403] [261 371]) - old, for 4 anim each
    [pout_expnor_ratio, out_stderrp_expnor, ~, out_stderrp_nor]= sj_ztestprop2([sum(out_last2d_Ncorrtrials_nor) sum(out_last2d_Ntrials_nor)] , [sum(out_last2d_Ncorrtrials_exp) sum(out_last2d_Ntrials_exp)]);
    pout_expnor_ratio; % 0.0014
    out_stderrp_expnor; % 0.0358
    
     % ([395 481] [321 403]) - old, for 4 anim each
    [pout_connor_ratio, out_stderrp_connor]= sj_ztestprop2([sum(out_last2d_Ncorrtrials_con) sum(out_last2d_Ntrials_con)] , [sum(out_last2d_Ncorrtrials_nor) sum(out_last2d_Ntrials_nor)]);
    pout_connor_ratio; % 0.0014
    out_stderrp_connor; % 0.0358
    
    errorbar(2, sum(out_last2d_Ncorrtrials_con)./sum(out_last2d_Ntrials_con),out_stderrp_con,'bd','MarkerSize',16,'Linewidth',3);
    errorbar(1, sum(out_last2d_Ncorrtrials_exp)./sum(out_last2d_Ntrials_exp),out_stderrp_exp,'ro','MarkerSize',16,'Linewidth',3);
    errorbar(3, sum(out_last2d_Ncorrtrials_nor)./sum(out_last2d_Ntrials_nor),out_stderrp_nor,'ks','MarkerSize',16,'Linewidth',3);
    if pout_expcon_ratio <0.05
        plot(2, sum(out_last2d_Ncorrtrials_con)./sum(out_last2d_Ntrials_con)+2*out_stderrp_con, 'r*','MarkerSize',16,'LineWidth',2);
    end
    if pout_expnor_ratio <0.05
        plot(3, sum(out_last2d_Ncorrtrials_nor)./sum(out_last2d_Ntrials_nor)+2*out_stderrp_nor, 'r*','MarkerSize',16,'LineWidth',2);
    end   
    axis([0 4 0 1.1]);
    set(gca,'YLim',[0.6 1])
    set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    title(['Outbound Performance in last ' num2str(lastdays) ' days: ztest ratio'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');

    % Save
    if savefig1==1,
        figfile = [figdir,'OutboundRawPerf_Last2Days_ztestratio'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 6.)Outbound Performance on extra days and pre-extra days - USE RAW PERF
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    %     % Raw Performance
    %     %a) Exper
    %     plot(0.9*ones(size(out_ex_perfday_expall)),expre_out_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',4);
    %     plot(1.1*ones(size(out_ex_perfday_expall)),ex_out_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',4)
    %     %b) Con
      plot(1.9*ones(size(out_ex_perfday_conall)),expre_out_rawperfday_con,'bd','MarkerSize',16,'Linewidth',3);
      plot(2.1*ones(size(out_ex_perfday_conall)),ex_out_rawperfday_con,'bd','MarkerSize',16,'Linewidth',3);
      for i=1:length(out_ex_perfday_conall)
          plot([1.9 2.1], [expre_out_rawperfday_con(i) ex_out_rawperfday_con(i)],'b-','Linewidth',2);
      end
    
    % Fitted Performance
    %a) Exper
%     plot(0.9*ones(size(out_ex_perfday_expall)),expre_out_perfday_exp,'ro','MarkerSize',16,'Linewidth',3);
%     plot(1.1*ones(size(out_ex_perfday_expall)),ex_out_perfday_exp,'ro','MarkerSize',16,'Linewidth',3)
    %b) Con
%     plot(1.9*ones(size(out_ex_perfday_conall)),expre_out_perfday_con,'bd','MarkerSize',16,'Linewidth',3);
%     plot(2.1*ones(size(out_ex_perfday_conall)),ex_out_perfday_con,'bd','MarkerSize',16,'Linewidth',3);
    
    % Lines between markers
%     for i=1:length(out_ex_perfday_expall)
%         plot([0.9 1.1], [expre_out_perfday_exp(i) ex_out_perfday_exp(i)],'r-','Linewidth',3);
%     end
    % Lines between markers
%     for i=1:length(out_ex_perfday_conall)
%         plot([1.9 2.1], [expre_out_perfday_con(i) ex_out_perfday_con(i)],'b-','Linewidth',3);
%     end
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.6 0.9]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Disruption';'Control'});
    title(['Outbound Performance before and after disruption'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pout_con_beforeafter = ranksum(expre_out_perfday_con,ex_out_perfday_con);
    
    % Save
    if savefig1==1,
        figfile = [figdir,'OutboundRawPerf_ExtraDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 7.) Outbound Ratio on extra days and pre-extra days - CAN ONLY USE RAW PERF
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Statistics
    
    % ([262 327] [229 317])
    [pout_con_ratio, out_stderrp_conbefaft,out_stderrp_conbef,out_stderrp_conaft]= sj_ztestprop2([sum(expre_out_Ncorrtrials_con) sum(expre_out_Ntrials_con)] , [sum(ex_out_Ncorrtrials_con) sum(ex_out_Ntrials_con)]);
    pout_con_ratio;
    out_stderrp_conbefaft;
    % ([182 265] [224 304])
    %[p_exp_ratio, stderrp_exp ]= sj_ztestprop2([sum(expre_out_Ncorrtrials_exp) sum(expre_out_Ntrials_exp)] , [sum(ex_out_Ncorrtrials_exp) sum(ex_out_Ntrials_exp)]);   
    
    %a) Exper
%     bar(0.5, sum(expre_out_Ncorrtrials_exp)./sum(expre_out_Ntrials_exp),'EdgeColor','r','Linewidth',3,'FaceColor','none');
%     errorbar(0.5, sum(expre_out_Ncorrtrials_exp)./sum(expre_out_Ntrials_exp),stderrp_exp,'r','Linewidth',3);
%     bar(1.5, sum(ex_out_Ncorrtrials_exp)./sum(ex_out_Ntrials_exp),'EdgeColor','r','Linewidth',3,'FaceColor','none');
%     errorbar(1.5, sum(ex_out_Ncorrtrials_exp)./sum(ex_out_Ntrials_exp),stderrp_exp,'r','Linewidth',3);
%     if p_exp_ratio <0.05
%         plot(1.5, sum(ex_out_Ncorrtrials_exp)./sum(ex_out_Ntrials_exp)+2*stderrp_exp, 'r*','MarkerSize',12,'LineWidth',2);
%     end
    %b) Con
    %bar(3.5, sum(expre_out_Ncorrtrials_con)./sum(expre_out_Ntrials_con),'EdgeColor','b','Linewidth',3,'FaceColor','none');
    errorbar(3.5, sum(expre_out_Ncorrtrials_con)./sum(expre_out_Ntrials_con),out_stderrp_conbef,'bd','MarkerSize',16,'Linewidth',3);
    %bar(4.5, sum(ex_out_Ncorrtrials_con)./sum(ex_out_Ntrials_con),'EdgeColor','b','Linewidth',3,'FaceColor','none');
    errorbar(4.5, sum(ex_out_Ncorrtrials_con)./sum(ex_out_Ntrials_con),out_stderrp_conaft,'bd','MarkerSize',16,'Linewidth',3);
    if pout_con_ratio <0.05
        plot(4.5, sum(ex_out_Ncorrtrials_con)./sum(ex_out_Ntrials_con)+2*out_stderrp_conbef, 'r*','MarkerSize',16,'LineWidth',2);
    end
        
    %axis([-0.5 5.5 0 1.05]);
    axis([2.8 5.3 0 1]);
    set(gca,'YLim',[0.6 0.91])
    set(gca,'XTick',[1 4],'XTickLabel',{'Disruption';'Control'});
    title(['Outbound Perf (ratio) for controls before and after disruption'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'OutboundPerfRawRatio_ExtraDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end % end if figout


%*******************************************************************
%*******************************************************************





if figin ==1
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1.) Inbound: Learning Trial, Session and Day and Ntrials
    % Trial
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    lt(1,:)=[mean(in_exp_lt),mean(in_con_lt)];
    %bar(lt);
    plot(2,in_con_lt(1),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,in_con_lt(2),'bd','MarkerSize',16,'Linewidth',3);
    plot(1.9,in_con_lt(3),'bd','MarkerSize',16,'Linewidth',3);
    %plot(2.1,in_con_lt(4),'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(2.1,in_con_lt(4),'bd','MarkerSize',16,'Linewidth',3);
    
    plot(1,in_exp_lt(1),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,in_exp_lt(2),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,in_exp_lt(3),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,in_exp_lt(4),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,in_exp_lt(5),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,in_exp_lt(6),'ro','MarkerSize',16,'Linewidth',3);
    
    %errorbar(0.8,mean(in_con_lt),sem(in_con_lt),'bd','MarkerSize',16,'Linewidth',4);
    %errorbar(1.2,mean(in_exp_lt),sem(in_exp_lt),'rd','MarkerSize',16,'Linewidth',4)
    
    axis([0 3 -5 200]);
    set(gca,'XTick',[1 2]);
    title('Inbound Learning Trial','FontSize',tfont,'Fontweight','normal');
    ylabel('Trials to reach Learning Criterion','FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');

    if ~isempty(normidx)
        for i=1:length(in_nor_lt)
            jitter(i)=0+0.05*(i-1);
        end
        %plot(3*ones(size(in_nor_lt))+jitter,in_nor_lt,'ks','MarkerSize',16,'Linewidth',4,'MarkerFaceColor','w');
        plot(2.9,in_nor_lt(1),'ks','MarkerSize',16,'Linewidth',3);
        plot(3,in_nor_lt(2),'ks','MarkerSize',16,'Linewidth',3);
        plot(3.1,in_nor_lt(3),'ks','MarkerSize',16,'Linewidth',3);
        plot(3,in_nor_lt(4),'ks','MarkerSize',16,'Linewidth',3);
        axis([0 4 -5 205]);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    end
    
    % Stats
    pin_expcon_lt = ranksum(in_con_lt,in_exp_lt); % p=0.8857
    pin_expnor_lt = ranksum(in_nor_lt,in_exp_lt); % p=1
    pin_connor_lt = ranksum(in_con_lt,in_nor_lt); % p=0.57
    pin_expvsconnor_lt = ranksum(in_exp_lt,[in_con_lt,in_nor_lt]); %0.9164

    % Save
    if savefig1==1,
        figfile = [figdir,'0Inbound_LearningTrial'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 2.) Day/Sess
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    ls(1,:)=[mean(in_exp_ls),mean(in_con_ls)];
    ld(1,:)=[mean(in_exp_ld),mean(in_con_ld)];
    %bar(lt);
    plot(2,in_con_ld(1),'bd','MarkerSize',16,'Linewidth',3);
    plot(2.1,in_con_ld(2),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,in_con_ld(3),'bd','MarkerSize',16,'Linewidth',3);
    plot(2.1,in_con_ld(4),'bd','MarkerSize',16,'Linewidth',3);
    
    plot(1.1,in_exp_ld(1),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,in_exp_ld(2),'ro','MarkerSize',16,'Linewidth',3);
    plot(0.9,in_exp_ld(3),'ro','MarkerSize',16,'Linewidth',3);
    plot(0.8,in_exp_ld(4),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,in_exp_ld(5),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,in_exp_ld(6),'ro','MarkerSize',16,'Linewidth',3);
    %errorbar(0.8,mean(in_con_ld),sem(in_con_ld),'bd','MarkerSize',16,'Linewidth',2);
    %errorbar(1.2,mean(in_exp_ld),sem(in_exp_ld),'rd','MarkerSize',16,'Linewidth',2)
    
    axis([0 3 0.5 8.5]);
    set(gca,'XTick',[1 2]);
    title('Inbound Learning Day','FontSize',tfont,'Fontweight','normal');
    ylabel('Days to reach Learning Criterion','FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    if ~isempty(normidx)
        for i=1:length(in_nor_ld)
            jitter(i)=0+0.1*(i-1);
        end
        jitter=[-0.1 0 0 0];
        plot(3*ones(size(in_nor_ld))+jitter,in_nor_ld,'ks','MarkerSize',16,'Linewidth',3);
        axis([0 4 0.5 8.5]);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    end
    
    % Stats
    pin_expcon_ld = ranksum(in_con_ld,in_exp_ld); %p=0.67
    pin_expnor_ld = ranksum(in_nor_ld,in_exp_ld); %p=0.67
    pin_connor_ld = ranksum(in_con_ld,in_nor_ld); %p=1
    pin_expvsconnor_ld = ranksum(in_exp_ld,[in_con_ld,in_nor_ld]); %p=0.48

    % Save
    if savefig1==1,
        figfile = [figdir,'0Inbound_LearningDay'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 3.) Total Ntrials in day_cutoff days
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    in_TotalNtr_expall =sum(in_ntrajspd_expall,2);
    in_TotalNtr_conall =sum(in_ntrajspd_conall,2);
    in_TotalNtr_norall =sum(in_ntrajspd_norall,2);
    %plot(in_TotalNtr_conall,2*ones(size(in_TotalNtr_conall)),'bd','MarkerSize',16,'Linewidth',4,'MarkerFaceColor','b');
    plot(2,in_TotalNtr_conall(1),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,in_TotalNtr_conall(2),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,in_TotalNtr_conall(3),'bd','MarkerSize',16,'Linewidth',3);
    plot(2,in_TotalNtr_conall(4),'bd','MarkerSize',16,'Linewidth',3);
    %plot(in_TotalNtr_expall,1*ones(size(in_TotalNtr_expall)),'ro','MarkerSize',16,'Linewidth',4);
    plot(1.07,in_TotalNtr_expall(1),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,in_TotalNtr_expall(2),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,in_TotalNtr_expall(3),'ro','MarkerSize',16,'Linewidth',3);
    plot(1,in_TotalNtr_expall(4),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,in_TotalNtr_expall(5),'ro','MarkerSize',16,'Linewidth',3);
    plot(1.2,in_TotalNtr_expall(6),'ro','MarkerSize',16,'Linewidth',3);
    
    axis([0 3 0 550]);
    set(gca,'XTick',[1 2]);
    title(['Total Number of Inbound Trials'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Number of Trials'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    if ~isempty(normidx)
        for i=1:length(in_TotalNtr_norall)
            jitter(i)=0+0.05*(i-1);
        end
        plot(3*ones(size(in_TotalNtr_norall))+jitter',in_TotalNtr_norall,'ks','MarkerSize',16,'Linewidth',3);
        axis([0 4 0 560]);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    end
    
    % Stats
    pin_expcon_Ntrials = ranksum(in_TotalNtr_conall,in_TotalNtr_expall);%p=0.48
    pin_expnor_Ntrials = ranksum(in_TotalNtr_norall,in_TotalNtr_expall);%p=0.35
    pin_connor_Ntrials = ranksum(in_TotalNtr_conall,in_TotalNtr_norall);%p=0.8857
    pin_expvsconnor_Ntrials = ranksum(in_TotalNtr_expall,[in_TotalNtr_conall;in_TotalNtr_norall]); %0.28
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0Inbound_NTrials'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 4.) Inbound Ntrials vs days
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    %     plot(in_ntrajspd_exp,'r.-','Linewidth',2,'MarkerSize',18);
    %     jbfill([1:length(in_ntrajvsd_exp)],in_ntrajvsd_expuperr,in_ntrajvsd_explowerr,'r','r',1,0.3);
    %     plot(in_ntrajspd_con,'b.-','Linewidth',2,'MarkerSize',18);
    %     jbfill([1:length(in_ntrajvsd_con)],in_ntrajvsd_conuperr,in_ntrajvsd_conlowerr,'b','b',1,0.3);
    errorbar([1:8],in_ntrajspd_exp,in_ntrajspd_experr,'ro-','Linewidth',3,'MarkerSize',12);
    errorbar([1.1:8.1],in_ntrajspd_con,in_ntrajspd_conerr,'bd-','Linewidth',3,'MarkerSize',12);
    [panovarm_Ntrdays_ConExpin] = anova_rm({in_ntrajspd_expall in_ntrajspd_conall},'off');
    %ConExp: pgrp=0.2424, ptime=0, pint=0.6034
    [panovarm_Ntrdays_ConExpNorin] = anova_rm({in_ntrajspd_expall in_ntrajspd_conall in_ntrajspd_norall},'off');
    %ConExpNor: pgrp=0.3908, ptime=0, pint=0.2867
%     if ~isempty(normidx)
%         errorbar([0.9:7.9],in_ntrajspd_nor,in_ntrajspd_norerr,'k.-','Linewidth',2,'MarkerSize',20);
%     end
    
    % Make Plot presentable
    title('Number of Inbound Trials vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Days','FontSize',xfont,'Fontweight','normal');
    ylabel('Inbound Trials','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 8.5 0 80]);
    set(gca,'XTick',[1:8],'XTickLabel',{[1:8]'});
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'0ConvsExp_Inbound_NTrialsPerDayVsDay'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 5.) Inbound Performance in last 2 days - USE RAW PERF ONLY
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    %     To Use Raw Performance rather than fitted performance, use in_rawperfday_con
    %     Raw Performance
        plot(2,in_rawperfday_con(1),'bd','MarkerSize',16,'Linewidth',3);
        plot(2.1,in_rawperfday_con(2),'bd','MarkerSize',16,'Linewidth',3);
        plot(2,in_rawperfday_con(3),'bd','MarkerSize',16,'Linewidth',3);
        plot(1.9,in_rawperfday_con(4),'bd','MarkerSize',16,'Linewidth',3);
        for i=1:length(in_rawperfday_exp)
            jitter(i)=0+0.05*(i-1);
        end
        %plot(1*ones(size(in_rawperfday_exp))+jitter,in_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',3);
        plot(1,in_rawperfday_exp(1),'ro','MarkerSize',16,'Linewidth',3);
        plot(1,in_rawperfday_exp(2),'ro','MarkerSize',16,'Linewidth',3);
        plot(0.9,in_rawperfday_exp(3),'ro','MarkerSize',16,'Linewidth',3);
        plot(0.9,in_rawperfday_exp(4),'ro','MarkerSize',16,'Linewidth',3);
        plot(1.2,in_rawperfday_exp(5),'ro','MarkerSize',16,'Linewidth',3);
        plot(1.2,in_rawperfday_exp(6),'ro','MarkerSize',16,'Linewidth',3);
    
      % Fitted Performance
%     plot(2,in_perfday_con(1),'bd','MarkerSize',16,'Linewidth',3);
%     plot(2.1,in_perfday_con(2),'bd','MarkerSize',16,'Linewidth',3);
%     plot(2,in_perfday_con(3),'bd','MarkerSize',16,'Linewidth',3);
%     plot(1.9,in_perfday_con(4),'bd','MarkerSize',16,'Linewidth',3);
%     plot(1*ones(size(in_perfday_expall)),in_perfday_exp,'ro','MarkerSize',16,'Linewidth',3);
    
    axis([0 3 0 1.05]);
    set(gca,'XTick',[1 2]);
    title(['Inbound Performance in last ' num2str(lastdays) ' days'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    if ~isempty(normidx)
        jitter=[];
        for i=1:length(in_rawperfday_nor)
            jitter(i)=0+0.05*(i-1);
        end
        % Raw Performance
        plot(3*ones(size(in_rawperfday_nor))+jitter,in_perfday_nor,'ks','MarkerSize',16,'Linewidth',3);
        % Fitted Performance
%       plot(3*ones(size(in_perfday_nor))+jitter,in_perfday_nor,'ks','MarkerSize',16,'Linewidth',3);
        
        axis([0 4 0 1.05]);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    end
    
    %Stats
    pin_expcon_last2days = ranksum(in_rawperfday_con,in_rawperfday_exp); %p=0.48
    pin_expnor_last2days = ranksum(in_rawperfday_nor,in_rawperfday_exp); %p=0.91
    pin_connor_last2days = ranksum(in_rawperfday_nor,in_rawperfday_con); %p=1
    pin_connor_last2days = ranksum(in_rawperfday_exp, [in_rawperfday_nor,in_rawperfday_con]); %p=0.35

    % Save
    if savefig1==1,
        figfile = [figdir,'0InboundRawPerf_Last2Days'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % 5b.) Inbound Performance in last 2 days - Plot Sem instead of points
    % for all 3 groups - use Raw Perf Only
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % To Use Raw Performance, use in_rawperfday_con
    errorbar(2,mean(in_rawperfday_con), sem(in_rawperfday_con),'bd','MarkerSize',16,'Linewidth',3);
    errorbar(1,mean(in_rawperfday_exp), sem(in_rawperfday_exp),'ro','MarkerSize',16,'Linewidth',3);
    errorbar(3,mean(in_rawperfday_nor), sem(in_rawperfday_nor),'ks','MarkerSize',16,'Linewidth',3);
 
    axis([0 4 0 1.05]);
    set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    title(['Inbound Performance in last ' num2str(lastdays) ' days'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');

    % Save
    if savefig1==1,
        figfile = [figdir,'InboundRawPerf_Last2Days_witherr'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 5c.)  Prop Correct ztest Ratio: Inbound Performance in last 2 days - 
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
     % ([481 496] [372 392]) - old, for 4 anim each
     %
    [pin_expcon_ratio, in_stderrp_expcon,in_stderrp_exp,in_stderrp_con]= sj_ztestprop2([sum(in_last2d_Ncorrtrials_con) sum(in_last2d_Ntrials_con)] , [sum(in_last2d_Ncorrtrials_exp) sum(in_last2d_Ntrials_exp)]);
    pin_expcon_ratio; % 0.1141
    in_stderrp_expcon; % 0.0131
    
      % ([401 432] [372 392]) - old, for 4 anim each
      %
    [pin_expnor_ratio, in_stderrp_expnor,~,in_stderrp_nor]= sj_ztestprop2([sum(in_last2d_Ncorrtrials_nor) sum(in_last2d_Ntrials_nor)] , [sum(in_last2d_Ncorrtrials_exp) sum(in_last2d_Ntrials_exp)]);
    pin_expnor_ratio; % 0.2173
    in_stderrp_expnor; % 0.0168
    
      % ([401 432] [372 392]) - old, for 4 anim each
      %
    [pin_connor_ratio, in_stderrp_connor]= sj_ztestprop2([sum(in_last2d_Ncorrtrials_con) sum(in_last2d_Ntrials_con)] , [sum(in_last2d_Ncorrtrials_nor) sum(in_last2d_Ntrials_nor)]);
    pin_connor_ratio; % 0.067
    in_stderrp_connor; % 0.0143
    
    errorbar(2, sum(in_last2d_Ncorrtrials_con)./sum(in_last2d_Ntrials_con),in_stderrp_con,'bd','MarkerSize',16,'Linewidth',3);
    errorbar(1, sum(in_last2d_Ncorrtrials_exp)./sum(in_last2d_Ntrials_exp),in_stderrp_exp,'ro','MarkerSize',16,'Linewidth',3);
    errorbar(3, sum(in_last2d_Ncorrtrials_nor)./sum(in_last2d_Ntrials_nor),in_stderrp_nor,'ks','MarkerSize',16,'Linewidth',3);
    if pin_expcon_ratio <0.05
        plot(2, sum(in_last2d_Ncorrtrials_con)./sum(in_last2d_Ntrials_con)+2*in_stderrp_expcon, 'r*','MarkerSize',16,'LineWidth',2);
    end
    if pin_expnor_ratio <0.05
        plot(3, sum(in_last2d_Ncorrtrials_nor)./sum(in_last2d_Ntrials_nor)+2*in_stderrp_expnor, 'r*','MarkerSize',16,'LineWidth',2);
    end   
    axis([0 4 0 1.1]);
    set(gca,'YLim',[0.6 1])
    set(gca,'XTick',[1 2 3],'XTickLabel',{'Disruption';'Control';'Normal'});
    title(['Inbound Performance in last ' num2str(lastdays) ' days: ztest ratio'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');

    % Save
    if savefig1==1,
        figfile = [figdir,'InboundRawPerf_Last2Days_ztestratio'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % 6.) Inbound Performance on extra days and pre-extra days - USE RAW PERF
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    %     % Raw Performance
    %     %a) Exper
    %     plot(0.9*ones(size(in_ex_perfday_expall)),expre_in_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',4);
    %     plot(1.1*ones(size(in_ex_perfday_expall)),ex_in_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',4)
    %     %b) Con
      plot(1.9*ones(size(in_ex_perfday_conall)),expre_in_rawperfday_con,'bd','MarkerSize',16,'Linewidth',3);
      plot(2.1*ones(size(in_ex_perfday_conall)),ex_in_rawperfday_con,'bd','MarkerSize',16,'Linewidth',3);
      for i=1:length(in_ex_perfday_conall)
          plot([1.9 2.1], [expre_in_rawperfday_con(i) ex_in_rawperfday_con(i)],'b-','Linewidth',2);
      end
    
    % Fitted Performance
    %a) Exper
%     plot(0.9*ones(size(in_ex_perfday_expall)),expre_in_perfday_exp,'ro','MarkerSize',16,'Linewidth',3);
%     plot(1.1*ones(size(in_ex_perfday_expall)),ex_in_perfday_exp,'ro','MarkerSize',16,'Linewidth',3)
    %b) Con
%     plot(1.9*ones(size(in_ex_perfday_conall)),expre_in_perfday_con,'bd','MarkerSize',16,'Linewidth',3);
%     plot(2.1*ones(size(in_ex_perfday_conall)),ex_in_perfday_con,'bd','MarkerSize',16,'Linewidth',3);
    
    % Lines between markers
%     for i=1:length(in_ex_perfday_expall)
%         plot([0.9 1.1], [expre_in_perfday_exp(i) ex_in_perfday_exp(i)],'r-','Linewidth',3);
%     end
    % Lines between markers
%     for i=1:length(in_ex_perfday_conall)
%         plot([1.9 2.1], [expre_in_perfday_con(i) ex_in_perfday_con(i)],'b-','Linewidth',3);
%     end
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.8 1.1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Disruption';'Control'});
    title(['Inbound Performance before and after disruption'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pin_con_beforeafter = ranksum(expre_in_perfday_con,ex_in_perfday_con);
    
    % Save
    if savefig1==1,
        figfile = [figdir,'InboundRawPerf_ExtraDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 7.) Inbound Ratio on extra days and pre-extra days - CAN ONLY USE RAW PERF
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Statistics
    
    % ([262 317] [229 317])
    [pin_con_ratio, in_stderrp_conbefaft, in_stderrp_conbef, in_stderrp_conaft]= sj_ztestprop2([sum(expre_in_Ncorrtrials_con) sum(expre_in_Ntrials_con)] , [sum(ex_in_Ncorrtrials_con) sum(ex_in_Ntrials_con)]);
    pin_con_ratio;
    in_stderrp_conbefaft;
    % ([182 265] [224 304])
    %[p_exp_ratio, stderrp_exp ]= sj_ztestprop2([sum(expre_in_Ncorrtrials_exp) sum(expre_in_Ntrials_exp)] , [sum(ex_in_Ncorrtrials_exp) sum(ex_in_Ntrials_exp)]);   
    
    %a) Exper
%     bar(0.5, sum(expre_in_Ncorrtrials_exp)./sum(expre_in_Ntrials_exp),'EdgeColor','r','Linewidth',3,'FaceColor','none');
%     errorbar(0.5, sum(expre_in_Ncorrtrials_exp)./sum(expre_in_Ntrials_exp),stderrp_exp,'r','Linewidth',3);
%     bar(1.5, sum(ex_in_Ncorrtrials_exp)./sum(ex_in_Ntrials_exp),'EdgeColor','r','Linewidth',3,'FaceColor','none');
%     errorbar(1.5, sum(ex_in_Ncorrtrials_exp)./sum(ex_in_Ntrials_exp),stderrp_exp,'r','Linewidth',3);
%     if p_exp_ratio <0.05
%         plot(1.5, sum(ex_in_Ncorrtrials_exp)./sum(ex_in_Ntrials_exp)+2*stderrp_exp, 'r*','MarkerSize',12,'LineWidth',2);
%     end
    %b) Con
    %bar(3.5, sum(expre_in_Ncorrtrials_con)./sum(expre_in_Ntrials_con),'EdgeColor','b','Linewidth',3,'FaceColor','none');
    errorbar(3.5, sum(expre_in_Ncorrtrials_con)./sum(expre_in_Ntrials_con),in_stderrp_conbef,'bd','MarkerSize',16,'Linewidth',3);
    %bar(4.5, sum(ex_in_Ncorrtrials_con)./sum(ex_in_Ntrials_con),'EdgeColor','b','Linewidth',3,'FaceColor','none');
    errorbar(4.5, sum(ex_in_Ncorrtrials_con)./sum(ex_in_Ntrials_con),in_stderrp_conaft,'bd','MarkerSize',16,'Linewidth',3);
    if pin_con_ratio <0.05
        plot(4.5, sum(ex_in_Ncorrtrials_con)./sum(ex_in_Ntrials_con)+2*in_stderrp_conaft, 'r*','MarkerSize',16,'LineWidth',2);
    end
        
    %axis([-0.5 5.5 0 1.05]);
    axis([2.8 5.3 0 1]);
    set(gca,'YLim',[0.85 1.05])
    set(gca,'XTick',[1 4],'XTickLabel',{'Disruption';'Control'});
    title(['Inbound Perf (ratio) for controls before and after disruption'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    % Save
    if savefig1==1,
        figfile = [figdir,'InboundPerfRawRatio_ExtraDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
            
   
    
end   % end if figin




%*******************************************************************
%*******************************************************************





%% Save Data
if savedata1==1,
    
    %savefile = '/data25/sjadhav/RippleInterruption/ProcessedData/BehSumm1';
    %savefile = '/data25/sjadhav/RippleInterruption/Figures';
    %savefile = sprintf('%s/ProcessedData/%s_LTP_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
    savefile = 'RipDis_BehSumm.mat';
    save(savefile,'Ripdis_summ');
end

%/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Behavior
x=1;

keyboard;
ranksum(out_con_lt,out_exp_lt);
ranksum(out_con_ld,out_exp_ld);
ranksum(out_perfday_con,out_perfday_exp); ranksum(out_perfday_nor,out_perfday_exp);
ranksum(out_rawperfday_con,out_rawperfday_exp); ranksum(out_rawperfday_nor,out_rawperfday_exp);
ranksum(expre_out_perfday_exp,ex_out_perfday_exp);
ranksum(expre_out_perfday_con,ex_out_perfday_con);




