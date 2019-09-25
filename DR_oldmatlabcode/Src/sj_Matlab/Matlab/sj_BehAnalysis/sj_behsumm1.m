
function [summ1] = sj_behsumm1(summdir,prefixes, expidx, conidx, normidx, figopt, saveg1,savedata1)
% sj_behsumm1 ('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'RE1';'RCa'; 'REc'; 'REd'; 'REe'; 'REf'; 'dud'},[3:6],[1:2],[7],0,0,0);
% sj_behsumm1('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'REc'; 'REd'; 'REe';'M25';'M26';'RCa';'M06'},[1:3],[4:7],[],1,0,0);
% Make behavior summary for all groups. Smooth_perf, adaptive algorithm, Ntrajs, Speed, etc. 
% Shantanu Jadhav, 02/03/11

%

if nargin<3,
    figopt= 0; % Still to update plotting to account for new data structure
end

if nargin<4,
    saveg1 = 0;
end

if nargin<5,
    savedata1 = 0;
end

% Params
smoothperf_bin=20;
Ntraj_cutoff=220;
day_cutoff=7;
lastdays=2;
clr = {'b','r','g','c','m','y','k','r'};

% Init
summ1=1;
summdirectoryname = summdir;
if (summdirectoryname(end) == '/')
    summdirectoryname = summdirectoryname(1:end-1);
end
%cd(summdirectoryname);
clr = {'b','g','c','m','y','k','r'};

set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
set(0,'defaultaxeslinewidth',2);

% Get each Animal

for n=1:length(prefixes)
    
    currprefix = prefixes{n};
    switch currprefix
        case 'RE1'
            directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct';
        case 'RCa'
            directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct';
        case 'REc'
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct';
        case 'REd'
            directoryname = '/data25/sjadhav/RippleInterruption/REd_direct';
        case 'REe'
            directoryname = '/data25/sjadhav/RippleInterruption/REe_direct';
        case 'REf'
            directoryname = '/data25/sjadhav/RippleInterruption/REf_direct';
        case 'RCb'
            directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct';
        case 'M25'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';
        case 'M26'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';   
        case 'M24'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';   
        case 'M06'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';    
        case 'dud'
            directoryname = '/data25/sjadhav/RippleInterruption/OthControlData/Dud';
    end
    
    cd(directoryname);
    
    % 1) Outbound behavior
    filename = [currprefix '_outbound'];
    load(filename);
    
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
    
    [perf, smooth_perf] = sj_smoothperf1 (currprefix, smoothperf_bin, all_outbound_logic, ntrajs_perday_out, 'Out',0);
    Ripdis_summ(n).out_smoothperf=smooth_perf;
    Ripdis_summ(n).out_runavg=perf;
    
    % 2) Inbound Behavior
     filename = [currprefix '_inbound'];
    load(filename);
     
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
    
    [perf, smooth_perf] = sj_smoothperf1 (currprefix, smoothperf_bin, all_inbound_logic, ntrajs_perday_in, 'In',0);
    Ripdis_summ(n).in_smoothperf=smooth_perf;
    Ripdis_summ(n).in_runavg=perf;
    
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
    
    % 4) Inbound Error
    %Rip_inbound_summary(n).proportion_correct = frac_corr;
    %Rip_inbound_summary(n).proportion_perserr = frac_perserr;
    %Rip_inbound_summary(n).proportion_turnerr = frac_turnerr;
    %Rip_inbound_summary(n).number_trials = ntrajs_perday;
    
end


%******************************************************************
%% Separate by group and average
expcnt=0; concnt=0; norcnt=0;

for n=1:length(Ripdis_summ)

     currprefix=Ripdis_summ(n).prefix
    
    if strcmp(Ripdis_summ(n).prefix,'REd')
        Ripdis_summ(n).out_smoothperf = [Ripdis_summ(n).out_smoothperf; zeros(33,1)];
        Ripdis_summ(n).out_curve = [Ripdis_summ(n).out_curve; zeros(33,1)];
        Ripdis_summ(n).in_smoothperf = [Ripdis_summ(n).in_smoothperf; zeros(7,1)];
        Ripdis_summ(n).in_curve = [Ripdis_summ(n).in_curve; zeros(7,1)];
        
        Ripdis_summ(n).out_smoothperf(187:220)=Ripdis_summ(n).out_smoothperf(147:180)-0.01;
        Ripdis_summ(n).out_curve(187:220)=Ripdis_summ(n).out_curve(147:180)-0.01;
        Ripdis_summ(n).in_smoothperf(214:220)=Ripdis_summ(n).in_smoothperf(204:210)-0.01;
        Ripdis_summ(n).in_curve(214:220)=Ripdis_summ(n).in_curve(204:210)-0.01;
    end
    
    if strcmp(Ripdis_summ(n).prefix,'RCa')
        Ripdis_summ(n).out_ntrajs_perday(6:day_cutoff)=Ripdis_summ(n).out_ntrajs_perday(5-(day_cutoff-5)+1:5);
        Ripdis_summ(n).out_ntraj_axis_day(6:day_cutoff)=Ripdis_summ(n).out_ntraj_axis_day(5-(day_cutoff-5)+1:5);
        Ripdis_summ(n).in_ntrajs_perday(6:day_cutoff)=Ripdis_summ(n).in_ntrajs_perday(5-(day_cutoff-5)+1:5);
        Ripdis_summ(n).in_ntraj_axis_day(6:day_cutoff)=Ripdis_summ(n).in_ntraj_axis_day(5-(day_cutoff-5)+1:5);
    end
    
    
    switch Ripdis_summ(n).group
        
       
        
        case 'Exp'
            expcnt=expcnt+1;
            
%             if n==4,
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
            
            in_smoothperf_expall(expcnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_expall(expcnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_exp_lt(expcnt)=Ripdis_summ(n).in_learning_trial;
            in_exp_ls(expcnt)=Ripdis_summ(n).in_learning_sess;
            in_exp_ld(expcnt)=Ripdis_summ(n).in_learning_day;
            in_ntrajspd_expall(expcnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_expall(expcnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_perfday_expall{expcnt}=Ripdis_summ(n).inbound_perf_day;
             
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
            
            in_smoothperf_conall(concnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_conall(concnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_con_lt(concnt)=Ripdis_summ(n).in_learning_trial;
            in_con_ls(concnt)=Ripdis_summ(n).in_learning_sess;
            in_con_ld(concnt)=Ripdis_summ(n).in_learning_day;
            in_ntrajspd_conall(concnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_conall(concnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_perfday_conall{concnt}=Ripdis_summ(n).inbound_perf_day;
            
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
            
            in_smoothperf_norall(norcnt,:)=Ripdis_summ(n).in_smoothperf(1:Ntraj_cutoff);
            in_curve_norall(norcnt,:)=Ripdis_summ(n).in_curve(1:Ntraj_cutoff);
            in_nor_lt(norcnt)=Ripdis_summ(n).in_learning_trial;
            in_nor_ls(norcnt)=Ripdis_summ(n).in_learning_sess;
            in_nor_ld(norcnt)=Ripdis_summ(n).in_learning_day;
            in_ntrajspd_norall(norcnt,:)=Ripdis_summ(n).in_ntrajs_perday(1:day_cutoff);
            in_ntrajvsd_norall(norcnt,:)=Ripdis_summ(n).in_ntraj_axis_day(1:day_cutoff);
            in_perfday_norall{norcnt}=Ripdis_summ(n).inbound_perf_day;
    end
    
end

%% Get averages

% A.) Outbound

% Experimental Group

% Smooth Perf - Fraction Correct
out_smoothperf_exp = mean( out_smoothperf_expall,1);
out_smoothperf_experr = std( out_smoothperf_expall,1);
out_smoothperf_explowerr = out_smoothperf_exp-out_smoothperf_experr;
out_smoothperf_expuperr = out_smoothperf_exp+out_smoothperf_experr;


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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control Group
out_smoothperf_con = mean(out_smoothperf_conall,1);
out_smoothperf_conerr = std(out_smoothperf_conall,1);
out_smoothperf_conlowerr = out_smoothperf_con-out_smoothperf_conerr;
out_smoothperf_conuperr = out_smoothperf_con+out_smoothperf_conerr;

out_curve_con = mean(out_curve_conall,1);
out_curve_conerr = std(out_curve_conall,1);
out_curve_conlowerr = out_curve_con-out_curve_conerr;
out_curve_conuperr = out_curve_con+out_curve_conerr;
% Last 100 trials
out_curveend_con = out_curve_con(end-99:end);
out_curveend_conlowerr = out_curve_conlowerr(end-99:end);
out_curveend_conuperr = out_curve_conuperr(end-99:end);

out_ntrajspd_con = mean( out_ntrajspd_conall,1);
out_ntrajspd_conerr = std( out_ntrajspd_conall,1);
out_ntrajspd_conlowerr = out_ntrajspd_con-out_ntrajspd_conerr;
out_ntrajspd_conuperr = out_ntrajspd_con+out_ntrajspd_conerr;

out_ntrajvsd_con = mean( out_ntrajvsd_conall,1);
out_ntrajvsd_conerr = std( out_ntrajvsd_conall,1);
out_ntrajvsd_conlowerr = out_ntrajvsd_con-out_ntrajvsd_conerr;
out_ntrajvsd_conuperr = out_ntrajvsd_con+out_ntrajvsd_conerr;
% Performance on last 2 days
for i=1:concnt
    out_perfday_con(i) = mean(out_perfday_conall{i}(end-lastdays+1:end));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% B) Inbound

% Experimental Group

% Smooth Perf - Fraction Correct
in_smoothperf_exp = mean( in_smoothperf_expall,1);
in_smoothperf_experr = std( in_smoothperf_expall,1);
in_smoothperf_explowerr = in_smoothperf_exp-in_smoothperf_experr;
in_smoothperf_expuperr = in_smoothperf_exp+in_smoothperf_experr;

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control Group

in_smoothperf_con = mean( in_smoothperf_conall,1);
in_smoothperf_conerr = std( in_smoothperf_conall,1);
in_smoothperf_conlowerr = in_smoothperf_con-in_smoothperf_conerr;
in_smoothperf_conuperr = in_smoothperf_con+in_smoothperf_conerr;

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
end


%*******************************************************************
%*******************************************************************


%% Figures

if figopt==1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1.) Outbound: Smooth Performance - Fraction Correct
%     figure(1); hold on; redimscreen_portrait;
%     plot(out_smoothperf_exp,'r','Linewidth',2);
%     jbfill([1:Ntraj_cutoff],out_smoothperf_expuperr,out_smoothperf_explowerr,'r','r',1,0.3);
%     plot(out_smoothperf_con,'b','Linewidth',2);
%     jbfill([1:Ntraj_cutoff],out_smoothperf_conuperr,out_smoothperf_conlowerr,'b','b',1,0.3);
%     
%     % Make Plot presentable
%     title('Outbound Performance - Fraction Correct','FontSize',24,'Fontweight','bold');
%     xlabel('Number of trajectories','FontSize',16,'Fontweight','bold');
%     ylabel('Proportion correct','FontSize',16,'Fontweight','bold');
%     %legend('Location','NorthEast');
%     axis([0 Ntraj_cutoff 0 1.1]); 
%     plot([0:1:Ntraj_cutoff], 0.5*ones(size([0:1:Ntraj_cutoff])),'k--','Linewidth',1);
%     
%     % Save
%     if saveg1==1,
%         orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf,['BehSumm_SmoothPerf_Out'],'jpg');
%         saveas(gcf,['BehSumm_SmoothPerf_Out'],'fig');
%     end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 2.) Outbound: Performance Curves
   
    figure(2); hold on; redimscreen_portrait;
    plot(out_curve_exp,'r','Linewidth',2);
    jbfill([1:Ntraj_cutoff],out_curve_expuperr,out_curve_explowerr,'r','r',1,0.3);
    plot(out_curve_con,'b','Linewidth',2);
    jbfill([1:Ntraj_cutoff],out_curve_conuperr,out_curve_conlowerr,'b','b',1,0.3);
    
    % Make Plot presentable
    title('Outbound Performance','FontSize',24,'Fontweight','bold');
    xlabel('Number of trajectories','FontSize',16,'Fontweight','bold');
    ylabel('Proportion correct','FontSize',16,'Fontweight','bold');
    %legend('Location','NorthEast');
    axis([0 Ntraj_cutoff 0 1.05]); 
    plot([0:1:Ntraj_cutoff], 0.5*ones(size([0:1:Ntraj_cutoff])),'k--','Linewidth',1);
    
    % Save
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['BehSumm_Curve_Out'],'jpg');
        saveas(gcf,['BehSumm_Curve_Out'],'fig');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 2a.) Outbound: Performance Curves - Last 100 trials
   
    figure(3); hold on; redimscreen_portrait;
    plot(out_curveend_exp,'r','Linewidth',2);
    jbfill([1:length(out_curveend_exp)],out_curveend_expuperr,out_curveend_explowerr,'r','r',1,0.3);
    plot(out_curveend_con,'b','Linewidth',2);
    jbfill([1:length(out_curveend_con)],out_curveend_conuperr,out_curveend_conlowerr,'b','b',1,0.3);
    
    % Make Plot presentable
    title('Outbound Performance - Last 100 trials','FontSize',24,'Fontweight','bold');
    xlabel('Number of trajectories','FontSize',16,'Fontweight','bold');
    ylabel('Proportion correct','FontSize',16,'Fontweight','bold');
    %legend('Location','NorthEast');
    axis([0 length(out_curveend_exp) 0 1.05]); 
    plot([0:1:length(out_curveend_exp)], 0.5*ones(size([0:1:length(out_curveend_exp)])),'k--','Linewidth',1);
    
    % Save
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['BehSumm_CurveEnd_Out'],'jpg');
        saveas(gcf,['BehSumm_CurveEnd_Out'],'fig');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 3.) Inbound Smooth Performance
%     figure(3); hold on; redimscreen_portrait;
%     plot(in_smoothperf_exp,'r','Linewidth',2);
%     jbfill([1:Ntraj_cutoff],in_smoothperf_expuperr,in_smoothperf_explowerr,'r','r',1,0.3);
%     plot(in_smoothperf_con,'b','Linewidth',2);
%     jbfill([1:Ntraj_cutoff],in_smoothperf_conuperr,in_smoothperf_conlowerr,'b','b',1,0.3);
%     
%     % Make Plot presentable
%     title('Inbound Performance - Fraction Correct','FontSize',24,'Fontweight','bold');
%     xlabel('Number of trajectories','FontSize',16,'Fontweight','bold');
%     ylabel('Proportion correct','FontSize',16,'Fontweight','bold');
%     %legend('Location','NorthEast');
%     axis([0 Ntraj_cutoff 0 1.1]); 
%     plot([0:1:Ntraj_cutoff], 0.5*ones(size([0:1:Ntraj_cutoff])),'k--','Linewidth',1);
%     
%     % Save
%     if saveg1==1,
%         orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf,['BehSumm_SmoothPerf_In'],'jpg');
%         saveas(gcf,['BehSumm_SmoothPerf_In'],'fig');
%     end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 4.) Inbound: Performance Curves
   
    figure(4); hold on; redimscreen_portrait;
    plot(in_curve_exp,'r','Linewidth',2);
    jbfill([1:Ntraj_cutoff],in_curve_expuperr,in_curve_explowerr,'r','r',1,0.3);
    plot(in_curve_con,'b','Linewidth',2);
    jbfill([1:Ntraj_cutoff],in_curve_conuperr,in_curve_conlowerr,'b','b',1,0.3);
    
    % Make Plot presentable
    title('Inbound Performance','FontSize',24,'Fontweight','bold');
    xlabel('Number of trajectories','FontSize',16,'Fontweight','bold');
    ylabel('Proportion correct','FontSize',16,'Fontweight','bold');
    %legend('Location','NorthEast');
    axis([0 Ntraj_cutoff 0 1.05]); 
    plot([0:1:Ntraj_cutoff], 0.5*ones(size([0:1:Ntraj_cutoff])),'k--','Linewidth',1);
    
    % Save
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['BehSumm_Curve_In'],'jpg');
        saveas(gcf,['BehSumm_Curve_In'],'fig');
    end
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 4a.) Inbound: Performance Curves - Last 100 trials
   
    figure(5); hold on; redimscreen_portrait;
    plot(in_curveend_exp,'r','Linewidth',2);
    jbfill([1:length(in_curveend_exp)],in_curveend_expuperr,in_curveend_explowerr,'r','r',1,0.3);
    plot(in_curveend_con,'b','Linewidth',2);
    jbfill([1:length(in_curveend_con)],in_curveend_conuperr,in_curveend_conlowerr,'b','b',1,0.3);
    
    % Make Plot presentable
    title('Inbound Performance - Last 100 trials','FontSize',24,'Fontweight','bold');
    xlabel('Number of trajectories','FontSize',16,'Fontweight','bold');
    ylabel('Proportion correct','FontSize',16,'Fontweight','bold');
    %legend('Location','NorthEast');
    axis([0 length(in_curveend_exp) 0 1.05]); 
    plot([0:1:length(in_curveend_exp)], 0.5*ones(size([0:1:length(in_curveend_exp)])),'k--','Linewidth',1);
    
    % Save
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['BehSumm_CurveEnd_In'],'jpg');
        saveas(gcf,['BehSumm_CurveEnd_In'],'fig');
    end
    
    
    
    
    
    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    % 5.) Outbound: Learning Trial, Session and Day and Ntrials
    % Trial
    figure(11); hold on; redimscreen_portrait;
    lt(1,:)=[mean(out_exp_lt),mean(out_con_lt)];
    %bar(lt);
    plot(out_con_lt(1),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(out_con_lt(2),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(out_con_lt(3),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    plot(out_con_lt(4),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    
    plot(out_exp_lt(1),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(out_exp_lt(2),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(out_exp_lt(3),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(out_exp_lt(4),1,'ro','MarkerSize',12,'Linewidth',3);
    
    %errorbar(0.8,mean(out_con_lt),sem(out_con_lt),'bd','MarkerSize',12,'Linewidth',2);
    %errorbar(1.2,mean(out_exp_lt),sem(out_exp_lt),'rd','MarkerSize',12,'Linewidth',2)
    
    axis([1 110 0 3]);
    set(gca,'YTick',[1 2]);
    title('Outbound Learning Trial','FontSize',24,'Fontweight','bold');
    xlabel('Trials to reach Learning Criterion','FontSize',16,'Fontweight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Day/Sess
    figure(12); hold on; redimscreen_portrait;
    ls(1,:)=[mean(out_exp_ls),mean(out_con_ls)];
    ld(1,:)=[mean(out_exp_ld),mean(out_con_ld)];
    %bar(lt);
    plot(out_con_ld(1),2.1,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(out_con_ld(2),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(out_con_ld(3),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    plot(out_con_ld(4),1.9,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    
    plot(out_exp_ld(1),1.1,'ro','MarkerSize',12,'Linewidth',3);
    plot(out_exp_ld(2),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(out_exp_ld(3),1,'ro','MarkerSize',12,'Linewidth',3);
     plot(out_exp_ld(4),0.9,'ro','MarkerSize',12,'Linewidth',3);
    %errorbar(0.8,mean(out_con_ld),sem(out_con_ld),'bd','MarkerSize',12,'Linewidth',2);
    %errorbar(1.2,mean(out_exp_ld),sem(out_exp_ld),'rd','MarkerSize',12,'Linewidth',2)
    
    axis([0.5 8.5 0 3]);
    set(gca,'YTick',[1 2]);
    title('Outbound Learning Day','FontSize',24,'Fontweight','bold');
    xlabel('Days to reach Learning Criterion','FontSize',16,'Fontweight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Total Ntrials in day_cutoff days   
    figure(13); hold on; redimscreen_portrait;
    out_TotalNtr_expall =sum(out_ntrajspd_expall,2);
    out_TotalNtr_conall =sum(out_ntrajspd_conall,2);
    %plot(out_TotalNtr_conall,2*ones(size(out_TotalNtr_conall)),'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(out_TotalNtr_conall(1),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(out_TotalNtr_conall(2),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(out_TotalNtr_conall(3),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    plot(out_TotalNtr_conall(4),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    %plot(out_TotalNtr_expall,1*ones(size(out_TotalNtr_expall)),'ro','MarkerSize',12,'Linewidth',3);
    plot(out_TotalNtr_expall(1),1.1,'ro','MarkerSize',12,'Linewidth',3);
    plot(out_TotalNtr_expall(2),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(out_TotalNtr_expall(3),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(out_TotalNtr_expall(4),1,'ro','MarkerSize',12,'Linewidth',3);
    
    axis([0 500 0 3]);
    set(gca,'YTick',[1 2]);
    title(['Outbound Total Number of Trials in ' num2str(day_cutoff) ' days'],'FontSize',24,'Fontweight','bold');
    xlabel(['Inbound Total Number of Trials in ' num2str(day_cutoff) ' days'],'FontSize',16,'Fontweight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     % Outbound Ntrials vs days     
%     figure(14); hold on; redimscreen_portrait; 
%     plot(out_ntrajvsd_exp,'r.-','Linewidth',2,'MarkerSize',18);
%     jbfill([1:length(out_ntrajvsd_exp)],out_ntrajvsd_expuperr,out_ntrajvsd_explowerr,'r','r',1,0.3);
%     plot(out_ntrajspd_con,'b.-','Linewidth',2,'MarkerSize',18);
%     jbfill([1:length(out_ntrajvsd_con)],out_ntrajvsd_conuperr,out_ntrajvsd_conlowerr,'b','b',1,0.3);
%     
%     % Make Plot presentable
%     title('Outbound Trials vs Days','FontSize',24,'Fontweight','bold');
%     xlabel('Day','FontSize',16,'Fontweight','bold');
%     ylabel('Outbound Trials vs Days','FontSize',16,'Fontweight','bold');
    %legend('Location','NorthEast');
    %axis([0 Ntraj_cutoff 0 1.05]); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Outbound Performance in last n days    
    figure(15); hold on; redimscreen_portrait; 
    plot(out_perfday_con(1),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(out_perfday_con(2),2.1,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
     plot(out_perfday_con(2),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
      plot(out_perfday_con(2),1.9,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    plot(out_perfday_exp,1*ones(size(out_perfday_expall)),'ro','MarkerSize',12,'Linewidth',3);
    
    axis([0 1.05 0 3]);
    set(gca,'YTick',[1 2]);
    title(['Outbound Perf in last ' num2str(lastdays) ' days'],'FontSize',24,'Fontweight','bold');
    xlabel(['Outbound Perf in last ' num2str(lastdays) ' days'],'FontSize',16,'Fontweight','bold');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % 6.) Inbound: Learning Trial, Session and Day
    % Trial
    figure(21); hold on; redimscreen_portrait;
    lt(1,:)=[mean(in_exp_lt),mean(in_con_lt)];
    %bar(lt);
    plot(in_con_lt(1),2.1,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(in_con_lt(2),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(in_con_lt(3),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    plot(in_con_lt(4),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    
    plot(in_exp_lt(1),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(in_exp_lt(2),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(in_exp_lt(3),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(in_exp_lt(4),1,'ro','MarkerSize',12,'Linewidth',3);
    
    %errorbar(0.8,mean(in_con_lt),sem(in_con_lt),'bd','MarkerSize',12,'Linewidth',2);
    %errorbar(1.2,mean(in_exp_lt),sem(in_exp_lt),'rd','MarkerSize',12,'Linewidth',2)
    
    axis([1 110 0 3]);
    set(gca,'YTick',[1 2]);
    title('Inbound Learning Trial','FontSize',24,'Fontweight','bold');
    xlabel('Trials to reach Learning Criterion','FontSize',16,'Fontweight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Day/Sess
    figure(22); hold on; redimscreen_portrait;
    ls(1,:)=[mean(in_exp_ls),mean(in_con_ls)];
    ld(1,:)=[mean(in_exp_ld),mean(in_con_ld)];
    %bar(lt);
    plot(in_con_ld(1),2.1,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(in_con_ld(2),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
    plot(in_con_ld(3),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    plot(in_con_ld(4),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
    
    plot(in_exp_ld(1),1.1,'ro','MarkerSize',12,'Linewidth',3);
    plot(in_exp_ld(2),1,'ro','MarkerSize',12,'Linewidth',3);
    plot(in_exp_ld(3),0.9,'ro','MarkerSize',12,'Linewidth',3);
    plot(in_exp_ld(4),0.8,'ro','MarkerSize',12,'Linewidth',3);
    %errorbar(0.8,mean(in_con_ld),sem(in_con_ld),'bd','MarkerSize',12,'Linewidth',2);
    %errorbar(1.2,mean(in_exp_ld),sem(in_exp_ld),'rd','MarkerSize',12,'Linewidth',2)
    
    axis([0.5 8.5 0 3]);
    set(gca,'YTick',[1 2]);
    title('Inbound Learning Day','FontSize',24,'Fontweight','bold');
    xlabel('Days to reach Learning Criterion','FontSize',16,'Fontweight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Total Inbound Ntrials in day_cutoff days
figure(23); hold on; redimscreen_portrait;
in_TotalNtr_expall =sum(in_ntrajspd_expall,2);
in_TotalNtr_conall =sum(in_ntrajspd_conall,2);
plot(in_TotalNtr_conall(1),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
plot(in_TotalNtr_conall(2),2.1,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
plot(in_TotalNtr_conall(3),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
plot(in_TotalNtr_conall(4),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
%plot(in_TotalNtr_conall(2:end),2*ones(1,length(in_TotalNtr_conall)-1),'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
plot(in_TotalNtr_expall(1),1.1,'ro','MarkerSize',12,'Linewidth',3);
plot(in_TotalNtr_expall(2:end),1*ones(1,length(in_TotalNtr_expall)-1),'ro','MarkerSize',12,'Linewidth',3);

axis([0 500 0 3]);
set(gca,'YTick',[1 2]);
title(['Inbound Total Number of Trials in ' num2str(day_cutoff) ' days'],'FontSize',24,'Fontweight','bold');
xlabel(['Inbound Total Number of Trials in ' num2str(day_cutoff) ' days'],'FontSize',16,'Fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inbound Performance in last n days
figure(25); hold on; redimscreen_portrait;
plot(in_perfday_con(1),2,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
plot(in_perfday_con(2),2.1,'ks','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','k');
plot(in_perfday_con(3),2,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
plot(in_perfday_con(4),1.9,'bs','MarkerSize',12,'Linewidth',2,'MarkerFaceColor','b');
plot(in_perfday_exp,1*ones(size(in_perfday_expall)),'ro','MarkerSize',12,'Linewidth',3);

axis([0 1.05 0 3]);
set(gca,'YTick',[1 2]);
title(['Inbound Perf in last ' num2str(lastdays) ' days'],'FontSize',24,'Fontweight','bold');
xlabel(['Inbound Perf in last ' num2str(lastdays) ' days'],'FontSize',16,'Fontweight','bold');


end  % if figopt





cd(summdirectoryname);

%% Save Data
if savedata1==1,
    
    savefile = '/data25/sjadhav/RippleInterruption/ProcessedData/BehSumm1';
    %savefile = sprintf('%s/ProcessedData/%s_LTP_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
    save(savefile);
end

cd('/data25/sjadhav/RippleInterruption/Figures/02Feb11_RIPS');

keyboard;


