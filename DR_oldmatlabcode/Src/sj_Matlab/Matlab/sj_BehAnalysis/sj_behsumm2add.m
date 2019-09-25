
function [summ1] = sj_behsumm2add(prefixes, expidx, conidx, normidx, figoptcurve, figoptpts, savefig1,savefig2, savedata1)

% sj_behsumm2add({'REc'; 'REd'; 'REe'; 'REf'; 'REg'; 'REh'; 'RCa';'RCb';'RCc';'RCd';'RNa';'RNb';'RNc';'RNd'},[1:6],[7:10],[11:14],0,1,0,0,0);

% Additional figures for beh summ - mainly to do with extra days

% f=1; % Choice For Plotting All Performance Curves
% f2=0; % Choice For Plotting ot Omitting Smooth Performance Curves
% f3=0; % Choice for plotting Group Comparison Stats

% Shantanu Jadhav, 9 Sep 11
%

if nargin<5,
    figoptcurve = 0; % Choice for plotting Fitted Curves
end

if nargin<6,
    figoptpts = 0; % Choice for plotting Fitted Curves
end

if nargin<7,
    savefig1 = 0; % Save figures or not
end

if nargin<8,
    savefig2 = 0; % Save figures or not
end

if nargin<9,
    savedata1 = 0; % Save Data or not
end


savefig3=0; % Use manually to save stuff you want

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Behavior/ExtraDays/';
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
exdays = 2;

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
    
    currprefix = prefixes{n}
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
    Ripdis_summ(n).out_ntraj_axis_day=ntraj_axis_day;
    
    [~, smooth_perf, rawperfvsday] = sj_smoothperf1 (currprefix, smoothperf_bin, all_outbound_logic, ntrajs_perday_out, 'Out',0);
    % Note: perf = smooth_perf. No difference
    
    Ripdis_summ(n).out_smoothperf=smooth_perf;
    Ripdis_summ(n).out_rawperfvsday=rawperfvsday;
    
    if n<=nexp+ncon
        if ~( strcmp(currprefix,'RCa') | strcmp(currprefix,'REc') | strcmp(currprefix,'REg') | strcmp(currprefix,'REh'))
            Ripdis_summ(n).ex_outbound_rawperf_day=ex_outbound_perf_day;
            Ripdis_summ(n).ex_out_curve=ex_outbound_curve;
            Ripdis_summ(n).ex_outbound_logic=ex_outbound_logic;
            
            [~, smooth_perf, ex_rawperf_day] = sj_smoothperf1 (currprefix, smoothperf_bin, ex_outbound_logic, ex_ntrajs_perday, 'Out',0);
            Ripdis_summ(n).ex_out_smoothperf=smooth_perf;
            
            %Ripdis_summ(n).ex_out_curve_witherr=ex_outbound_curve_witherr;
            
        else
            Ripdis_summ(n).ex_outbound_rawperf_day=[];
            Ripdis_summ(n).ex_out_curve=[];
            Ripdis_summ(n).ex_outbound_logic=[];
            Ripdis_summ(n).ex_out_smoothperf=[];

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
    Ripdis_summ(n).in_ntraj_axis_day=ntraj_axis_day;
    Ripdis_summ(n).in_ntraj_axis_sess=ntraj_axis_sess;
    %Ripdis_summ(n).inbound_perf_day=inbound_perf_day;
    
    [~, smooth_perf, rawperfvsday] = sj_smoothperf1 (currprefix, smoothperf_bin, all_inbound_logic, ntrajs_perday_in, 'In',0);
    % Note: perf = smooth_perf. No difference
    
    Ripdis_summ(n).in_smoothperf=smooth_perf;
    Ripdis_summ(n).in_rawperfvsday=rawperfvsday;
    
    if n<=nexp+ncon
        if ~( strcmp(currprefix,'RCa') | strcmp(currprefix,'REc') | strcmp(currprefix,'REg') | strcmp(currprefix,'REh') ) 
            % No extra days for these animals
            Ripdis_summ(n).ex_inbound_rawperf_day=ex_inbound_perf_day;
            Ripdis_summ(n).ex_in_curve=ex_inbound_curve;
            Ripdis_summ(n).ex_inbound_logic=ex_inbound_logic;
            
            [~, smooth_perf, ex_rawperf_day] = sj_smoothperf1 (currprefix, smoothperf_bin, ex_inbound_logic, ex_ntrajs_perday, 'Out',0);
            Ripdis_summ(n).ex_in_smoothperf=smooth_perf;
            
            %Ripdis_summ(n).ex_in_curve_witherr=ex_inbound_curve_witherr;
            
        else
            Ripdis_summ(n).ex_inbound_rawperf_day=[];
            Ripdis_summ(n).ex_in_curve=[];
            Ripdis_summ(n).ex_inbound_logic=[];
            Ripdis_summ(n).ex_out_smoothperf=[];
            
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
            out_ntrajspd_expall(expcnt,:)=Ripdis_summ(n).out_ntrajs_perday(1:day_cutoff);
            out_ntrajvsd_expall(expcnt,:)=Ripdis_summ(n).out_ntraj_axis_day(1:day_cutoff);
            out_rawperfvsday_expall{expcnt}=Ripdis_summ(n).out_rawperfvsday;
            
            in_smoothperf_expall(expcnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_expall(expcnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_ntrajspd_expall(expcnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_expall(expcnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_rawperfvsday_expall{expcnt}=Ripdis_summ(n).in_rawperfvsday;
            
            if ~isempty(Ripdis_summ(n).ex_outbound_rawperf_day)
                expcnt_ex=expcnt_ex+1;
                
                % out extra 2 days
                out_ex_rawperfday_expall(expcnt_ex,:)=Ripdis_summ(n).ex_outbound_rawperf_day; 
                out_ex_curve_expall{expcnt_ex}=Ripdis_summ(n).ex_out_curve;
                out_ex_smoothperf_expall{expcnt_ex}=Ripdis_summ(n).ex_out_smoothperf;
                out_ex_logic_expall{expcnt_ex}=Ripdis_summ(n).ex_outbound_logic;
                
                % out last 2 days
                out_exbef_rawperfday_expall(expcnt_ex,:) = Ripdis_summ(n).out_rawperfvsday(end-exdays+1:end); 
                out_exbef_curve_expall{expcnt_ex}=Ripdis_summ(n).out_curve...
                    ( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: end);
                out_exbef_smoothperf_expall{expcnt_ex}=Ripdis_summ(n).out_smoothperf...
                    ( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: end);
                out_exbef_logic_expall{expcnt_ex}=Ripdis_summ(n).outbound_logic...
                    ( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: end);
                
                % in extra 2 days
                in_ex_rawperfday_expall(expcnt_ex,:)=Ripdis_summ(n).ex_inbound_rawperf_day; 
                in_ex_curve_expall{expcnt_ex}=Ripdis_summ(n).ex_in_curve;
                in_ex_smoothperf_expall{expcnt_ex}=Ripdis_summ(n).ex_in_smoothperf;
                in_ex_logic_expall{expcnt_ex}=Ripdis_summ(n).ex_inbound_logic;
                
                % in last 2 days
                in_exbef_rawperfday_expall(expcnt_ex,:) = Ripdis_summ(n).in_rawperfvsday(end-exdays+1:end); 
                in_exbef_curve_expall{expcnt_ex}=Ripdis_summ(n).in_curve...
                    ( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: end);
                in_exbef_smoothperf_expall{expcnt_ex}=Ripdis_summ(n).in_smoothperf...
                    ( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: end);
                in_exbef_logic_expall{expcnt_ex}=Ripdis_summ(n).inbound_logic...
                    ( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: end);
                
            end
            
        case 'Con'
            concnt=concnt+1;
            
            out_smoothperf_conall(concnt,:)=Ripdis_summ(n).out_smoothperf(1:Ntraj_cutoff);
            out_curve_conall(concnt,:)=Ripdis_summ(n).out_curve(1:Ntraj_cutoff);            
            out_ntrajspd_conall(concnt,:)=Ripdis_summ(n).out_ntrajs_perday(1:day_cutoff);
            out_ntrajvsd_conall(concnt,:)=Ripdis_summ(n).out_ntraj_axis_day(1:day_cutoff);
            out_rawperfvsday_conall{concnt}=Ripdis_summ(n).out_rawperfvsday;
            
            in_smoothperf_conall(concnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_conall(concnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_ntrajspd_conall(concnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_conall(concnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_rawperfvsday_conall{concnt}=Ripdis_summ(n).in_rawperfvsday;
            
            if ~isempty(Ripdis_summ(n).ex_outbound_rawperf_day)
                concnt_ex=concnt_ex+1;
                
                % out extra 2 days
                out_ex_rawperfday_conall(concnt_ex,:)=Ripdis_summ(n).ex_outbound_rawperf_day; 
                out_ex_curve_conall{concnt_ex}=Ripdis_summ(n).ex_out_curve;
                out_ex_smoothperf_conall{concnt_ex}=Ripdis_summ(n).ex_out_smoothperf;
                out_ex_logic_conall{concnt_ex}=Ripdis_summ(n).ex_outbound_logic;
                
                % out last 2 days
                out_exbef_rawperfday_conall(concnt_ex,:) = Ripdis_summ(n).out_rawperfvsday(end-exdays+1:end); 
                out_exbef_curve_conall{concnt_ex}=Ripdis_summ(n).out_curve...
                    ( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: end);
                out_exbef_smoothperf_conall{concnt_ex}=Ripdis_summ(n).out_smoothperf...
                    ( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: end);
                out_exbef_logic_conall{concnt_ex}=Ripdis_summ(n).outbound_logic...
                    ( Ripdis_summ(n).out_ntraj_axis_day(end-2)+1: end);
                
                % in extra 2 days
                in_ex_rawperfday_conall(concnt_ex,:)=Ripdis_summ(n).ex_inbound_rawperf_day; %extra 2 days
                in_ex_curve_conall{concnt_ex}=Ripdis_summ(n).ex_in_curve;
                in_ex_smoothperf_conall{concnt_ex}=Ripdis_summ(n).ex_in_smoothperf;
                in_ex_logic_conall{concnt_ex}=Ripdis_summ(n).ex_inbound_logic;
                
                % in last 2 days
                in_exbef_rawperfday_conall(concnt_ex,:) = Ripdis_summ(n).in_rawperfvsday(end-exdays+1:end); % last 2 days
                in_exbef_curve_conall{concnt_ex}=Ripdis_summ(n).in_curve...
                    ( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: end);
                in_exbef_smoothperf_conall{concnt_ex}=Ripdis_summ(n).in_smoothperf...
                    ( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: end);
                in_exbef_logic_conall{concnt_ex}=Ripdis_summ(n).inbound_logic...
                    ( Ripdis_summ(n).in_ntraj_axis_day(end-2)+1: end);
                
            end
            
        case 'Nor'
            norcnt=norcnt+1;
            
            out_smoothperf_norall(norcnt,:)=Ripdis_summ(n).out_smoothperf(1:Ntraj_cutoff);
            out_curve_norall(norcnt,:)=Ripdis_summ(n).out_curve(1:Ntraj_cutoff);
            out_ntrajspd_norall(norcnt,:)=Ripdis_summ(n).out_ntrajs_perday(1:day_cutoff);
            out_ntrajvsd_norall(norcnt,:)=Ripdis_summ(n).out_ntraj_axis_day(1:day_cutoff);
            out_rawperfvsday_norall{norcnt}=Ripdis_summ(n).out_rawperfvsday;
            
            in_smoothperf_norall(norcnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_norall(norcnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_ntrajspd_norall(norcnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_norall(norcnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_rawperfvsday_norall{norcnt}=Ripdis_summ(n).in_rawperfvsday;
    end
    
end

%% Get averages

%**************************************************
% A.) Outbound
%**************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experimental Group

for i=1:expcnt_ex
    
    % Raw perf
    ex_out_rawperfday_exp(i) = mean(out_ex_rawperfday_expall(i,:));
    exbef_out_rawperfday_exp(i) = mean(out_exbef_rawperfday_expall(i,:));
   
    % Curve: before after
    
    ex_out_curve_exp(i,:) = out_ex_curve_expall{i}(1:78); % No of trials are are 78, 130, 96
    if i==1, 
        exbef_out_curve_exp(i,:) = out_exbef_curve_expall{i}(1:61); 
    else
        exbef_out_curve_exp(i,:) = out_exbef_curve_expall{i}(end-61+1:end); % No of trials are are 71 (last 10 are 0), 102, 105
    end
    
    % Smoothperf: before after
    ex_out_smoothperf_exp(i,:) = out_ex_smoothperf_expall{i}(1:73); % No of trials are are 73, 125, 91
    if i==1, 
        exbef_out_smoothperf_exp(i,:) = out_exbef_smoothperf_expall{i}(1:61); 
    else
        exbef_out_smoothperf_exp(i,:) = out_exbef_smoothperf_expall{i}(end-61+1:end); % No of trials are are 66(last 5 are 0), 97, 100
    end
    
    % Ratio Test: Number of correct trials and no of total trials for 2 extra and 2 pre-extra days
    ex_out_Ntrials_exp(i) = length(out_ex_logic_expall{i});
    ex_out_Ncorrtrials_exp(i) = length(find(out_ex_logic_expall{i}==1));
    exbef_out_Ntrials_exp(i) = length(out_exbef_logic_expall{i});
    exbef_out_Ncorrtrials_exp(i) = length(find(out_exbef_logic_expall{i}==1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control Group


% Performance on extra days and Pre-extra days Limited to animals with
% Extra Days
for i=1:concnt_ex
    
    % Raw perf
    ex_out_rawperfday_con(i) = mean(out_ex_rawperfday_conall(i,:));
    exbef_out_rawperfday_con(i) = mean(out_exbef_rawperfday_conall(i,:));
   
    % Curve: before after
    
    ex_out_curve_con(i,:) = out_ex_curve_conall{i}(1:90); % No of trials are are 98, 129, 90
    exbef_out_curve_con(i,:) = out_exbef_curve_conall{i}(end-100+1:end); % No of trials are are 105, 122, 100
    
    % Smoothperf: before after
    ex_out_smoothperf_con(i,:) = out_ex_smoothperf_conall{i}(1:85); % No of trials are are 93, 124, 85
    exbef_out_smoothperf_con(i,:) = out_exbef_smoothperf_conall{i}(end-95+1:end); % No of trials are are 100, 117, 95
    
    % Ratio Test: Number of correct trials and no of total trials for 2 extra and 2 pre-extra days
    ex_out_Ntrials_con(i) = length(out_ex_logic_conall{i});
    ex_out_Ncorrtrials_con(i) = length(find(out_ex_logic_conall{i}==1));
    exbef_out_Ntrials_con(i) = length(out_exbef_logic_conall{i});
    exbef_out_Ncorrtrials_con(i) = length(find(out_exbef_logic_conall{i}==1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%**************************************************
% B) Inbound
%**************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experimental Group

for i=1:expcnt_ex
   
    % Raw perf
    ex_in_rawperfday_exp(i) = mean(in_ex_rawperfday_expall(i,:));
    exbef_in_rawperfday_exp(i) = mean(in_exbef_rawperfday_expall(i,:));
   
    % Curve: before after
    ex_in_curve_exp(i,:) = in_ex_curve_expall{i}(1:82); % No of trials are are 82, 137, 98
    % 3rd one has problems  - the 0s should be 1s
    if i==3, ex_in_curve_exp(i,38)=1;  ex_in_curve_exp(i,65)=1; end
    exbef_in_curve_exp(i,:) = in_exbef_curve_expall{i}(end-60+1:end); % No of trials are are 60, 113, 109
    
    % Smoothperf: before after
    ex_in_smoothperf_exp(i,:) = in_ex_smoothperf_expall{i}(1:77); % No of trials are are 77, 132, 93
    exbef_in_smoothperf_exp(i,:) = in_exbef_smoothperf_expall{i}(end-55+1:end); % No of trials are are 55, 108, 104
    
    % Ratio Test: Number of correct trials and no of total trials for 2 extra and 2 pre-extra days
    ex_in_Ntrials_exp(i) = length(in_ex_logic_expall{i});
    ex_in_Ncorrtrials_exp(i) = length(find(in_ex_logic_expall{i}==1));
    exbef_in_Ntrials_exp(i) = length(in_exbef_logic_expall{i});
    exbef_in_Ncorrtrials_exp(i) = length(find(in_exbef_logic_expall{i}==1));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control Group

% Performance on extra days and Pre-extra days Limited to animals with
% Extra Days
for i=1:concnt_ex
   
    % Raw perf
    ex_in_rawperfday_con(i) = mean(in_ex_rawperfday_conall(i,:));
    exbef_in_rawperfday_con(i) = mean(in_exbef_rawperfday_conall(i,:));
   
    % Curve: before after
    ex_in_curve_con(i,:) = in_ex_curve_conall{i}(1:85); % No of trials are are 108, 132, 85
    exbef_in_curve_con(i,:) = in_exbef_curve_conall{i}(end-100+1:end); % No of trials are are 114, 128, 100
    
    % Smoothperf: before after
    ex_in_smoothperf_con(i,:) = in_ex_smoothperf_conall{i}(1:80); % No of trials are are 103, 127, 80
    exbef_in_smoothperf_con(i,:) = in_exbef_smoothperf_conall{i}(end-72+1:end); % No of trials are are 72, 109, 85
    
    % Ratio Test: Number of correct trials and no of total trials for 2 extra and 2 pre-extra days
    ex_in_Ntrials_con(i) = length(in_ex_logic_conall{i});
    ex_in_Ncorrtrials_con(i) = length(find(in_ex_logic_conall{i}==1));
    exbef_in_Ntrials_con(i) = length(in_exbef_logic_conall{i});
    exbef_in_Ncorrtrials_con(i) = length(find(in_exbef_logic_conall{i}==1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%*******************************************************************
%*******************************************************************
%                          FIGURES
%*******************************************************************
%*******************************************************************

% Go to Figure Save Directory
cd(summdir);

if figoptcurve==1  % Choice For Plotting Fitted Curves
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1.) Outbound: Fitted Performance Curves before and after - Con
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end

    plot(mean(exbef_out_curve_con),'b','Linewidth',3);
    jbfill([1:length(exbef_out_curve_con)],mean(exbef_out_curve_con)+std(exbef_out_curve_con),...
        mean(exbef_out_curve_con)-std(exbef_out_curve_con),'b','b',1,1);
    
    out_exmean = mean(ex_out_curve_con);
    out_exerr = std(ex_out_curve_con);

    plot(length(exbef_out_curve_con)+1:length(exbef_out_curve_con)+length(ex_out_curve_con),...
        out_exmean,'c','Linewidth',3);
    jbfill(length(exbef_out_curve_con)+1:length(exbef_out_curve_con)+length(ex_out_curve_con),...
        out_exmean+out_exerr,...
        out_exmean-out_exerr,'c','c',1,1);
 
    plot(length(exbef_out_curve_con)*ones(size([0:0.1:1.1])),[0:0.1:1.1],'k--','Linewidth',1);

    % Make Plot presentable
    title('Outbound Curve Con: Bef-Aft','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    set(gca,'YLim',[0 1]);
    %axis([0 Ntraj_cutoff 0 1.01]);
    set(gca,'TickDir','out');

    % Save
    if savefig1==1,
        figfile = [figdir,'Con_BefAft_OutCurve'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 2.) Outbound: Smooth Performance Curves before and after - Con
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end

    plot(mean(exbef_out_smoothperf_con),'b','Linewidth',3);
    jbfill([1:length(exbef_out_smoothperf_con)],mean(exbef_out_smoothperf_con)+std(exbef_out_smoothperf_con),...
        mean(exbef_out_smoothperf_con)-std(exbef_out_smoothperf_con),'b','b',1,1);
    
    out_exmean = mean(ex_out_smoothperf_con);
    out_exerr = std(ex_out_smoothperf_con);
    
    plot(length(exbef_out_smoothperf_con)+1:length(exbef_out_smoothperf_con)+length(ex_out_smoothperf_con),...
        out_exmean,'c','Linewidth',3);
    jbfill(length(exbef_out_smoothperf_con)+1:length(exbef_out_smoothperf_con)+length(ex_out_smoothperf_con),...
        out_exmean+out_exerr,...
        out_exmean-out_exerr,'c','c',1,1);
 
    plot(length(exbef_out_smoothperf_con)*ones(size([0:0.1:1.1])),[0:0.1:1.1],'k--','Linewidth',1);

    % Make Plot presentable
    title('Outbound SmPerf Con: Bef-Aft','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    set(gca,'YLim',[0 1]);
    %axis([0 Ntraj_cutoff 0 1.01]);
    set(gca,'TickDir','out');

    % Save
%     if savefig1==1,
%         figfile = [figdir,'Con_BefAft_OutSmPerf'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    % 3.) Inbound: Fitted Performance Curves before and after - Con
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end

    plot(mean(exbef_in_curve_con),'b','Linewidth',3);
    jbfill([1:length(exbef_in_curve_con)],mean(exbef_in_curve_con)+std(exbef_in_curve_con),...
        mean(exbef_in_curve_con)-std(exbef_in_curve_con),'b','b',1,1);
    
    in_exmean = mean(ex_in_curve_con);
    in_exerr = std(ex_in_curve_con);
    
    plot(length(exbef_in_curve_con)+1:length(exbef_in_curve_con)+length(ex_in_curve_con),...
        in_exmean,'c','Linewidth',3);
    jbfill(length(exbef_in_curve_con)+1:length(exbef_in_curve_con)+length(ex_in_curve_con),...
        in_exmean+in_exerr,...
        in_exmean-in_exerr,'c','c',1,1);
 
    plot(length(exbef_in_curve_con)*ones(size([0:0.1:1.1])),[0:0.1:1.1],'k--','Linewidth',1);

    % Make Plot presentable
    title('Inbound Curve Con: Bef-Aft','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    set(gca,'YLim',[0 1.05]);
    %axis([0 Ntraj_cutoff 0 1.01]);
    set(gca,'TickDir','out');

    % Save
%     if savefig1==1,
%         figfile = [figdir,'Con_BefAft_InCurve'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 4.) Inbound: Smooth Performance Curves before and after - Con
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end

    plot(mean(exbef_in_smoothperf_con),'b','Linewidth',3);
    jbfill([1:length(exbef_in_smoothperf_con)],mean(exbef_in_smoothperf_con)+std(exbef_in_smoothperf_con),...
        mean(exbef_in_smoothperf_con)-std(exbef_in_smoothperf_con),'b','b',1,1);
    
    in_exmean = mean(ex_in_smoothperf_con);
    in_exerr = std(ex_in_smoothperf_con);
    
    plot(length(exbef_in_smoothperf_con)+1:length(exbef_in_smoothperf_con)+length(ex_in_smoothperf_con),...
        in_exmean,'c','Linewidth',3);
    jbfill(length(exbef_in_smoothperf_con)+1:length(exbef_in_smoothperf_con)+length(ex_in_smoothperf_con),...
        in_exmean+in_exerr,...
        in_exmean-in_exerr,'c','c',1,1);
 
    plot(length(exbef_in_smoothperf_con)*ones(size([0:0.1:1.1])),[0:0.1:1.1],'k--','Linewidth',1);

    % Make Plot presentable
    title('Inbound SmPerf Con: Bef-Aft','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    set(gca,'YLim',[0 1.05]);
    %axis([0 Ntraj_cutoff 0 1.01]);
    set(gca,'TickDir','out');

    % Save
    if savefig1==1,
        figfile = [figdir,'Con_BefAft_OutSmPerf'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
         
    % 5.) Outbound: Fitted Performance Curves before and after - Exp
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end

    plot(mean(exbef_out_curve_exp),'r','Linewidth',3);
    jbfill([1:length(exbef_out_curve_exp)],mean(exbef_out_curve_exp)+std(exbef_out_curve_exp),...
        mean(exbef_out_curve_exp)-std(exbef_out_curve_exp),'r','r',1,1);
    
    out_exmean = mean(ex_out_curve_exp);
    out_exerr = std(ex_out_curve_exp);

    plot(length(exbef_out_curve_exp)+1:length(exbef_out_curve_exp)+length(ex_out_curve_exp),...
        out_exmean,'m','Linewidth',3);
    jbfill(length(exbef_out_curve_exp)+1:length(exbef_out_curve_exp)+length(ex_out_curve_exp),...
        out_exmean+out_exerr,...
        out_exmean-out_exerr,'m','m',1,1);
 
    plot(length(exbef_out_curve_exp)*ones(size([0:0.1:1.1])),[0:0.1:1.1],'k--','Linewidth',1);

    % Make Plot presentable
    title('Outbound Curve Exp: Bef-Aft','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    set(gca,'YLim',[0 1]);
    %axis([0 Ntraj_cutoff 0 1.01]);
    set(gca,'TickDir','out');

    % Save
%     if savefig1==1,
%         figfile = [figdir,'Exp_BefAft_OutCurve'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 6.) Inbound: Fitted Performance Curves before and after - Exp
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end

    plot(mean(exbef_in_curve_exp),'r','Linewidth',3);
    jbfill([1:length(exbef_in_curve_exp)],mean(exbef_in_curve_exp)+std(exbef_in_curve_exp),...
        mean(exbef_in_curve_exp)-std(exbef_in_curve_exp),'r','r',1,1);
    
    in_exmean = mean(ex_in_curve_exp);
    in_exerr = std(ex_in_curve_exp);

    plot(length(exbef_in_curve_exp)+1:length(exbef_in_curve_exp)+length(ex_in_curve_exp),...
        in_exmean,'m','Linewidth',3);
    jbfill(length(exbef_in_curve_exp)+1:length(exbef_in_curve_exp)+length(ex_in_curve_exp),...
        in_exmean+in_exerr,...
        in_exmean-in_exerr,'m','m',1,1);
 
    plot(length(exbef_in_curve_exp)*ones(size([0:0.1:1.1])),[0:0.1:1.1],'k--','Linewidth',1);

    % Make Plot presentable
    title('Inbound Curve Exp: Bef-Aft','FontSize',tfont,'Fontweight','normal');
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal');
    ylabel('Probability Correct','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    set(gca,'YLim',[0 1.05]);
    %axis([0 Ntraj_cutoff 0 1.01]);
    set(gca,'TickDir','out');

    % Save
%     if savefig1==1,
%         figfile = [figdir,'Exp_BefAft_InCurve'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%*******************************************************************
%*******************************************************************








if figoptpts==1
        
    
    % After second review - 03Apr2012. Plot days 7,8 and 9,10 separately
    % a.) Con Outbound Performance on extra days and pre-extra days
    % With Individual Points and Lines between markers
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Plot Con
    plot([1,2,3,4],[out_exbef_rawperfday_conall(1,:),out_ex_rawperfday_conall(1,:)],'bd--','MarkerSize',16,'Linewidth',3);
    plot([1,2,3,4],[out_exbef_rawperfday_conall(2,:),out_ex_rawperfday_conall(2,:)],'bd:','MarkerSize',16,'Linewidth',3);
    plot([1,2,3,4],[out_exbef_rawperfday_conall(3,:),out_ex_rawperfday_conall(3,:)],'bd-','MarkerSize',16,'Linewidth',3);
    
    % Get all Exp Anim Last 2 days
    for i=1:length(out_rawperfvsday_expall)
        out_exp_last2d(i,:) = out_rawperfvsday_expall{i}(7:8);
    end
    for i=1:length(out_rawperfvsday_expall)
        plot([6,7],[out_exp_last2d(i,:)],'ro-','MarkerSize',16,'Linewidth',3);
    end
    
    axis([0 8 0.4 1.0]);
    set(gca,'YLim',[0.4 1]);
    set(gca,'XTick',[1 2 3 4 6 7],'XTickLabel',{'7','8','9','10','7','8'});
    %set(gca,'XTick',[1 2],'XTickLabel',{'Control'});
    title(['Con: OutPerf befaft disruption - Lines:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pout_con_beforeafter = ranksum(exbef_out_rawperfday_con,ex_out_rawperfday_con); % p=0.1
    
    % Save
    if savefig2==1,
        figfile = [figdir,'0Con_BefAft_OutRawPerf_Lines_axis2'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
       % b.) Con Inbound Performance on extra days and pre-extra days
    % With Individual Points and Lines between markers
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Plot Con
    plot([1,2,3,4],[in_exbef_rawperfday_conall(1,:),in_ex_rawperfday_conall(1,:)],'bd--','MarkerSize',16,'Linewidth',3);
    plot([1,2,3,4],[in_exbef_rawperfday_conall(2,:),in_ex_rawperfday_conall(2,:)],'bd:','MarkerSize',16,'Linewidth',3);
    plot([1,2,3,4],[in_exbef_rawperfday_conall(3,:),in_ex_rawperfday_conall(3,:)],'bd-','MarkerSize',16,'Linewidth',3);
    
    % Get all Exp Anim Last 2 days
    for i=1:length(in_rawperfvsday_expall)
        in_exp_last2d(i,:) = in_rawperfvsday_expall{i}(7:8);
    end
    for i=1:length(in_rawperfvsday_expall)
        plot([6,7],[in_exp_last2d(i,:)],'ro-','MarkerSize',16,'Linewidth',3);
    end
    
    axis([0 8 0.4 1.0]);
    set(gca,'YLim',[0.4 1]);
    set(gca,'XTick',[1 2 3 4 6 7],'XTickLabel',{'7','8','9','10','7','8'});
    %set(gca,'XTick',[1 2],'XTickLabel',{'Control'});
    title(['Con: OutPerf befaft disruption - Lines:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pin_con_beforeafter = ranksum(exbef_in_rawperfday_con,ex_in_rawperfday_con); % p=1
    
    % Save
    if savefig2==1,
        figfile = [figdir,'0Con_BefAft_InRawPerf_Lines_axis2'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
   % ------------------------------------------------
    % 1.) Con Outbound Performance on extra days and pre-extra days - USE RAW PERF
    % Same as in sj_behsumm2
    
    % With Individual Points and Lines between markers
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    plot(1.9*ones(size(exbef_out_rawperfday_con)),exbef_out_rawperfday_con,'bd','MarkerSize',16,'Linewidth',3);
    plot(2.1*ones(size(ex_out_rawperfday_con)),ex_out_rawperfday_con,'bd','MarkerSize',16,'Linewidth',3);
    for i=1:length(ex_out_rawperfday_con)
        plot([1.9 2.1], [exbef_out_rawperfday_con(i) ex_out_rawperfday_con(i)],'b-','Linewidth',2);
    end
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.5 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Control'});
    title(['Con: OutPerf befaft disruption - Lines:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pout_con_beforeafter = ranksum(exbef_out_rawperfday_con,ex_out_rawperfday_con); % p=0.1
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Con_BefAft_OutRawPerf_Lines_axis2'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 2.) With Mean and error before and after
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    errorbar(1.9,mean(exbef_out_rawperfday_con),sem(exbef_out_rawperfday_con),'bd','MarkerSize',16,'Linewidth',3);
    errorbar(2.1,mean(ex_out_rawperfday_con),sem(ex_out_rawperfday_con),'bd','MarkerSize',16,'Linewidth',3);
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.6 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Control'});
    title(['Con: OutPerf befaft disruption - Sem:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pout_con_beforeafter = ranksum(exbef_out_rawperfday_con,ex_out_rawperfday_con); % p=0.1
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Con_BefAft_OutRawPerf_WithErr'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 3.) Con - Outbound Ratio on extra days and pre-extra days - CAN ONLY USE RAW PERF
    % Same as sj_behsumm2
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Calc and Statistics
    
    % ([262 327] [229 317])
    [pout_con_ratio, out_stderrp_conbefaft, out_stderrp_conbef, out_stderrp_conaft]= sj_ztestprop2([sum(exbef_out_Ncorrtrials_con) sum(exbef_out_Ntrials_con)] , [sum(ex_out_Ncorrtrials_con) sum(ex_out_Ntrials_con)]);
    pout_con_ratio; % p=0.0188
    out_stderrp_conbefaft; % 0.0335
    %z=2.3498;
  
    errorbar(1.9, sum(exbef_out_Ncorrtrials_con)./sum(exbef_out_Ntrials_con),out_stderrp_conbef,'bd','MarkerSize',16,'Linewidth',3);
    errorbar(2.1, sum(ex_out_Ncorrtrials_con)./sum(ex_out_Ntrials_con),out_stderrp_conaft,'bd','MarkerSize',16,'Linewidth',3);
    if pout_con_ratio <0.05
        plot(2.1, sum(ex_out_Ncorrtrials_con)./sum(ex_out_Ntrials_con)+2*out_stderrp_conaft, 'r*','MarkerSize',16,'LineWidth',2);
    end
        
    %axis([-0.5 5.5 0 1.05]);
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.6 1])
    set(gca,'XTick',[2],'XTickLabel',{'Control'});
    title(['Outbound Perf (ratio) for controls before and after disruption'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Con_BefAft_OutPropRatio'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % 4.) Con Inbound Performance on extra days and pre-extra days - USE RAW PERF
    % Same as in sj_behsumm2
    
    % With Individual Points and Lines between markers
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    plot(1.9*ones(size(exbef_in_rawperfday_con)),exbef_in_rawperfday_con,'bd','MarkerSize',16,'Linewidth',3);
    plot(2.1*ones(size(ex_in_rawperfday_con)),ex_in_rawperfday_con,'bd','MarkerSize',16,'Linewidth',3);
    for i=1:length(ex_in_rawperfday_con)
        plot([1.9 2.1], [exbef_in_rawperfday_con(i) ex_in_rawperfday_con(i)],'b-','Linewidth',2);
    end
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.5 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Control'});
    title(['Con: InPerf befaft disruption - Lines:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pin_con_beforeafter = ranksum(exbef_in_rawperfday_con,ex_in_rawperfday_con); % p=0.1
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Con_BefAft_InRawPerf_Lines_axis2'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 5.) With Mean and error before and after
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    errorbar(1.9,mean(exbef_in_rawperfday_con),sem(exbef_in_rawperfday_con),'bd','MarkerSize',16,'Linewidth',3);
    errorbar(2.1,mean(ex_in_rawperfday_con),sem(ex_in_rawperfday_con),'bd','MarkerSize',16,'Linewidth',3);
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.6 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Control'});
    title(['Con: InPerf befaft disruption - Sem:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pin_con_beforeafter = ranksum(exbef_in_rawperfday_con,ex_in_rawperfday_con); % p=0.1
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Con_BefAft_InRawPerf_WithErr'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 6.) Con - Inbound Ratio on extra days and pre-extra days - CAN ONLY USE RAW PERF
    % Same as sj_behsumm2
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Calc and Statistics
    
    % ([328 342] [312 325])
    [pin_con_ratio, in_stderrp_conbefaft, in_stderrp_conbef, in_stderrp_conaft]= sj_ztestprop2([sum(exbef_in_Ncorrtrials_con) sum(exbef_in_Ntrials_con)] , [sum(ex_in_Ncorrtrials_con) sum(ex_in_Ntrials_con)]);
    pin_con_ratio; % p=0.9511
    in_stderrp_conbefaft; % 0.0153
    % z = 0.0613
  
    errorbar(1.9, sum(exbef_in_Ncorrtrials_con)./sum(exbef_in_Ntrials_con),in_stderrp_conbef,'bd','MarkerSize',16,'Linewidth',3);
    errorbar(2.1, sum(ex_in_Ncorrtrials_con)./sum(ex_in_Ntrials_con),in_stderrp_conaft,'bd','MarkerSize',16,'Linewidth',3);
    if pin_con_ratio <0.05
        plot(2.1, sum(ex_in_Ncorrtrials_con)./sum(ex_in_Ntrials_con)+2*in_stderrp_conaft, 'r*','MarkerSize',16,'LineWidth',2);
    end
        
    %axis([-0.5 5.5 0 1.05]);
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.6 1])
    set(gca,'XTick',[2],'XTickLabel',{'Control'});
    title(['Inbound Perf (ratio) for controls before and after disruption'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Con_BefAft_InPropRatio'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    % SIX FIGURES FOR EXPTAL
    
     % 1.) Exp Outbound Performance on extra days and pre-extra days - USE RAW PERF
    % Same as in sj_behsumm2
    
    % With Individual Points and Lines between markers
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    plot(1.9*ones(size(exbef_out_rawperfday_exp)),exbef_out_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',3);
    plot(2.1*ones(size(ex_out_rawperfday_exp)),ex_out_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',3);
    for i=1:length(ex_out_rawperfday_exp)
        plot([1.9 2.1], [exbef_out_rawperfday_exp(i) ex_out_rawperfday_exp(i)],'r-','Linewidth',2);
    end
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.5 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp'});
    title(['Exp: OutPerf befaft disruption - Lines:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pout_exp_beforeafter = ranksum(exbef_out_rawperfday_exp,ex_out_rawperfday_exp); % p=0.4
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Exp_BefAft_OutRawPerf_Lines'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 2.) With Mean and error before and after
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    errorbar(1.9,mean(exbef_out_rawperfday_exp),sem(exbef_out_rawperfday_exp),'ro','MarkerSize',16,'Linewidth',3);
    errorbar(2.1,mean(ex_out_rawperfday_exp),sem(ex_out_rawperfday_exp),'ro','MarkerSize',16,'Linewidth',3);
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.5 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp'});
    title(['Exp: OutPerf befaft disruption - Sem:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pout_exp_beforeafter = ranksum(exbef_out_rawperfday_exp,ex_out_rawperfday_exp); % p=0.4
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Exp_BefAft_OutRawPerf_WithErr'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 3.) Exp - Outbound Ratio on extra days and pre-extra days - CAN ONLY USE RAW PERF
    % Same as sj_behsumm2
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Calc and Statistics
    
    % ([182 265] [224 304])
    [pout_exp_ratio, out_stderrp_expbefaft, out_stderrp_expbef, out_stderrp_expaft]= sj_ztestprop2([sum(exbef_out_Ncorrtrials_exp) sum(exbef_out_Ntrials_exp)] , [sum(ex_out_Ncorrtrials_exp) sum(ex_out_Ntrials_exp)]);
    pout_exp_ratio; % p=0.1878
    out_stderrp_expbefaft; % 0.0380
    %z=1.3172
  
    errorbar(1.9, sum(exbef_out_Ncorrtrials_exp)./sum(exbef_out_Ntrials_exp),out_stderrp_expbef,'ro','MarkerSize',16,'Linewidth',3);
    errorbar(2.1, sum(ex_out_Ncorrtrials_exp)./sum(ex_out_Ntrials_exp),out_stderrp_expaft,'ro','MarkerSize',16,'Linewidth',3);
    if pout_exp_ratio <0.05
        plot(2.1, sum(ex_out_Ncorrtrials_exp)./sum(ex_out_Ntrials_exp)+2*out_stderrp_expaft, 'r*','MarkerSize',16,'LineWidth',2);
    end
        
    %axis([-0.5 5.5 0 1.05]);
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.5 1])
    set(gca,'XTick',[2],'XTickLabel',{'Exp'});
    title(['Outbound Perf (ratio) for exp before and after disruption'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Exp_BefAft_OutPropRatio'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % 4.) Exp Inbound Performance on extra days and pre-extra days - USE RAW PERF
    % Same as in sj_behsumm2
    
    % With Individual Points and Lines between markers
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    plot(1.9*ones(size(exbef_in_rawperfday_exp)),exbef_in_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',3);
    plot(2.1*ones(size(ex_in_rawperfday_exp)),ex_in_rawperfday_exp,'ro','MarkerSize',16,'Linewidth',3);
    for i=1:length(ex_in_rawperfday_exp)
        plot([1.9 2.1], [exbef_in_rawperfday_exp(i) ex_in_rawperfday_exp(i)],'r-','Linewidth',2);
    end
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.6 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp'});
    title(['Exp: InPerf befaft disruption - Lines:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pin_exp_beforeafter = ranksum(exbef_in_rawperfday_exp,ex_in_rawperfday_exp); % p=0.7
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Exp_BefAft_InRawPerf_Lines_axis2'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 5.) With Mean and error before and after
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    errorbar(1.9,mean(exbef_in_rawperfday_exp),sem(exbef_in_rawperfday_exp),'ro','MarkerSize',16,'Linewidth',3);
    errorbar(2.1,mean(ex_in_rawperfday_exp),sem(ex_in_rawperfday_exp),'ro','MarkerSize',16,'Linewidth',3);
    
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.5 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp'});
    title(['Exp: InPerf befaft disruption - Sem:All Anim'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Probability Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    %Stats
    pin_exp_beforeafter = ranksum(exbef_in_rawperfday_exp,ex_in_rawperfday_exp); % p=0.7
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Exp_BefAft_InRawPerf_WithErr'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end   
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 6.) Exp - Inbound Ratio on extra days and pre-extra days - CAN ONLY USE RAW PERF
    % Same as sj_behsumm2
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Calc and Statistics
    
    % ([264 282] [299 317])
    [pin_exp_ratio, in_stderrp_expbefaft, in_stderrp_expbef, in_stderrp_expaft]= sj_ztestprop2([sum(exbef_in_Ncorrtrials_exp) sum(exbef_in_Ntrials_exp)] , [sum(ex_in_Ncorrtrials_exp) sum(ex_in_Ntrials_exp)]);
    pin_exp_ratio; % p=0.7172
    in_stderrp_expbefaft; % 0.0195
    % z = 0.3622
  
    errorbar(1.9, sum(exbef_in_Ncorrtrials_exp)./sum(exbef_in_Ntrials_exp),in_stderrp_expbef,'ro','MarkerSize',16,'Linewidth',3);
    errorbar(2.1, sum(ex_in_Ncorrtrials_exp)./sum(ex_in_Ntrials_exp),in_stderrp_expaft,'ro','MarkerSize',16,'Linewidth',3);
    if pin_exp_ratio <0.05
        plot(2.1, sum(ex_in_Ncorrtrials_exp)./sum(ex_in_Ntrials_exp)+2*in_stderrp_expaft, 'r*','MarkerSize',16,'LineWidth',2);
    end
        
    %axis([-0.5 5.5 0 1.05]);
    axis([1.5 2.5 0 1.05]);
    set(gca,'YLim',[0.5 1])
    set(gca,'XTick',[2],'XTickLabel',{'Exp'});
    title(['Inbound Perf (ratio) for exp before and after disruption'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Proportion Correct'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'TickDir','out');
    
    % Save
    if savefig2==1,
        figfile = [figdir,'Exp_BefAft_InPropRatio'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        % ([182 265] [224 304])
    %[p_exp_ratio, stderrp_exp ]= sj_ztestprop2([sum(expre_out_Ncorrtrials_exp) sum(expre_out_Ntrials_exp)] , [sum(ex_out_Ncorrtrials_exp) sum(ex_out_Ntrials_exp)]);   


%*******************************************************************
%*******************************************************************


end % end figoptpts

    




   




%*******************************************************************
%*******************************************************************





%% Save Data
if savedata1==1,
    
    %savefile = '/data25/sjadhav/RippleInterruption/ProcessedData/BehSumm1';
    %savefile = '/data25/sjadhav/RippleInterruption/Figures';
    %savefile = sprintf('%s/ProcessedData/%s_LTP_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
    savefile = 'RipDis_BehSumm.mat';
    save(savefile,'Ripdis_summ_add');
end

%/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Behavior
x=1;

keyboard;





