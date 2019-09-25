
function [wellsdio, rewardinfo] = sj_findwellsfromdio2(directoryname,prefix,days,epochs, dionos, varargin)
%%% Shantanu - Fix errors

% [wellsdio] =
% sj_findwellsfromdio1('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',[1:3],[2 4],[26 27 28]);
% [wellsdio] = sj_findwellsfromdio1('/data25/sjadhav/RippleInterruption/REd_direct','REd',[1:2],[2 4],[22 23 24]);

% [wellsdio] = sj_findwellsfromdio2('/data25/sjadhav/HPExpt/HPb_direct','HPb',[1],[4 6],[24 23 22]);
% [wellsdio] = sj_findwellsfromdio2('/data25/sjadhav/HPExpt/HPb_direct','HPb',[2],[2 4],[24 23 22]);
% sj_findwellsfromdio1('/data25/sjadhav/HPExpt/HPa_direct','HPa',[2],[2 4],[24 23 22]);
%DR works% [wellsdio, rewardinfo] = sj_findwellsfromdio2('/mnt/data25/sjadhav/HPExpt/HPc', 'HPc',2,[2 4], [24 23 22]);
% [wellsdio, rewardinfo] = sj_findwellsfromdio2('/mnt/data25/sjadhav/HPExpt/HPc', 'HPc',2,[2 4], [24 23 22]);

%% From DIO, gets well start and end for all unique trajectories

lowercasethree = '';
%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
        
    end
end

for day=days,
    
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    % Check that center well has more triggers than outer wells by loading epoch2
    % Is this necessary? - NO - GET THE ORDER RIGHT
%     li=[1 2 3];
%     le1=length(DIO{day}{2}{dionos(1)}.pulsetimes);
%     le2=length(DIO{day}{2}{dionos(2)}.pulsetimes);
%     le3=length(DIO{day}{2}{dionos(3)}.pulsetimes);
%     [len,li] = sort([le1,le2,le3],'descend');
%     usedios=dionos(li);
    %li=[1 2 3];
    
    % IF you know order of DIOs - enter here. Center well first
    usedios = [27,26,25];
    outdios=[9,10,11];
    
    for epoch=epochs,
        
        %Init
        well_start=[]; well_end=[]; well_startend=[]; wellseq_curr=[]; welltrigtime_curr=[]; tmpwells=[];
        
        DIO1=DIO{day}{epoch}{usedios(1)}; nones = length(DIO1.pulsetimes); % Center well 1: Bits 8/24 in
        DIO2=DIO{day}{epoch}{usedios(2)}; nzeros = length(DIO2.pulsetimes); % Well 0: Bits 9/25 in
        DIO3=DIO{day}{epoch}{usedios(3)}; ntwos = length(DIO3.pulsetimes); % Well 2: Bits 10/26 in
        
        % Output wells
        DIOout1=DIO{day}{epoch}{outdios(1)}; nones_out = length(DIOout1.pulsetimes); % Center well 1: Bits 8 out/24
        DIOout2=DIO{day}{epoch}{outdios(2)}; nzeros_out = length(DIOout2.pulsetimes); % Well 0: Bits 9 out/25
        DIOout3=DIO{day}{epoch}{outdios(3)}; ntwos_out = length(DIOout3.pulsetimes); % Well 2: Bits 10 out/26
        out_trigt = [DIOout1.pulsetimes(:,1);DIOout2.pulsetimes(:,1);DIOout3.pulsetimes(:,1)];
        [out_trigtx, out_ord] = sort(out_trigt);
        out_wells = [1*ones(nones_out,1);0*ones(nzeros_out,1);2*ones(ntwos_out,1)];
        out_rewwells = out_wells(out_ord);
        
        
%         figure; hold on;
%         plot(DIO1.pulsetimes./10000,ones(size(DIO1.pulsetimes)),'bo');
%         plot(DIO2.pulsetimes./10000,0*ones(size(DIO2.pulsetimes)),'bo');
%         plot(DIO3.pulsetimes./10000,2*ones(size(DIO3.pulsetimes)),'bo');
%         
%         plot(DIOout1.pulsetimes./10000,ones(size(DIOout1.pulsetimes)),'rd');
%         plot(DIOout2.pulsetimes./10000,0*ones(size(DIOout2.pulsetimes)),'rd');
%         plot(DIOout3.pulsetimes./10000,2*ones(size(DIOout3.pulsetimes)),'rd');
        
        
        
        % At which wells did the trigger occur? Real trajectories are the ones
        % where the trigger well changed eg. 0 0 0 1  : transition from
        % trigger at well 0 to well 1
        % -------------------------------------------------------------
        wells = [1*ones(nones,1);0*ones(nzeros,1);2*ones(ntwos,1)];
        trigt = [DIO1.pulsetimes(:,1);DIO2.pulsetimes(:,1);DIO3.pulsetimes(:,1)];
        [tfx,ti] = sort(trigt); % sorted triggers
        wfx=wells(ti); % cosrted wells according to trigger times
        swit=find(diff(wfx)~=0)+1; % find where trajectories changed: start of trajectory
        well_start=[wfx(1);wfx(swit(1:end-1))]; % start well of trajectory - add back first well and get rid of last, since no trajectory for last well
        well_end=[wfx(swit)]; % end well of trajectory
        well_startend=[well_start,well_end]; % well start and end for detecting outbound and inbound
        % For saving
        wellsdio{day}{epoch}=well_startend; 
        
        
        
        % Now use the well sequence to get outbound and inbound as well as
        % rewarded or unrewarded
        % -------------------------------------------------------
        wellseq_curr=[wfx(1);wfx(swit)]; % Sequence of wells triggerred
        welltrigtime_curr=[tfx(1);tfx(swit)]; % Trigger times at each of these wells in sequence
        
        % IF First Well is 1, then first outbound was captured. Keep it aside temporarily
        % Remove first outbound if it exists. Put back later
        % Reason is sometimes first trigger at center well may not be
        % captured since data save began too late
        % -------------------------------------------
        tmpwells = well_startend; flag_firstout=0;
         if wellseq_curr(1)==1, 
             tmpwells(1,:)=[]; 
             flag_firstout=1;
         end
%         if wellseq_curr(1)~=1,
%             tmpseq=[1;tmpseq];
%             first = [1, well_startend(1,1)];
%             tmpwells = [first; tmpwells];
%         end
           
        % Use Sequence of wells to get outbound and inbound trajectories,
        % and correct/rewarded or incorrect/unrewarded
        % ---------------------------------------
        %%% Find outbound
        outbound_stidx = find(tmpwells(:,1)==1);
        outbound_wellstend = tmpwells(outbound_stidx,:);
        outbound_logic = zeros(length(outbound_wellstend),1);
        corridx = find( (tmpwells(outbound_stidx,2)~=tmpwells(outbound_stidx-1,1))...
            & (tmpwells(outbound_stidx,2)~=1) );
        correct_out = outbound_stidx(corridx);
        wrong_out = setdiff(outbound_stidx,correct_out);
        %%% Find inbound
        inbound_stidx = find(tmpwells(:,1)~=1);
        inbound_wellstend = tmpwells(inbound_stidx,:);
        % Find out which inbound is correct, and return the logic
        inbound_logic = zeros(length(inbound_wellstend),1);
        corridx = find(inbound_wellstend(:,2)==1);
        correct_in = inbound_stidx(corridx);
        wrong_in = setdiff(inbound_stidx,correct_in);
        %inbound_logic(corr) = 1;
        
        % IMP - Put first Outbound Back In If It Was Removed - Have to Push all indexes By 1
        % ------------------------------------------------------------------------
        %if wellseq_curr(1)~=1,
        if wellseq_curr(1)==1,
            first = well_startend(1,:);
            tmpwells = [first; tmpwells];           
            % Push out and in idxs by 1
            correct_out = correct_out+1; correct_out=[1;correct_out];
            correct_in = correct_in+1;
            wrong_out = wrong_out+1;
            wrong_in = wrong_in+1;
            outbound_stidx = outbound_stidx+1; outbound_stidx=[1;outbound_stidx]; 
            inbound_stidx = inbound_stidx+1;
        end
                
        % Combine the outbound and inbound together as well as correct or not
        % -----------------------------------------------------------
        all_correct = sort([correct_out; correct_in]);
        all_wrong = sort([wrong_out; wrong_in]);
        allwell_curr = sort([all_correct;all_wrong]);
        rewtime_curr = welltrigtime_curr(all_correct+1);  % Why +1. Because Reward at end of trajectory: next well in sequence
        norewtime_curr = welltrigtime_curr(all_wrong+1);
        rewardedwell_seq = well_startend(all_correct,2); unrewardedwell_seq = well_startend(all_wrong,2);
        wellseq_reconstructed = wellseq_curr(allwell_curr+1);
        
        % Now prepare the reward structure for outputting
        % ---------------------------------------------
        if wellseq_curr(1)==1,
            wellseq_out =  wellseq_curr(2:end); % Skipping 1st well (center), since that could have been captured in data acquisition or not
            welltrigtime_out = welltrigtime_curr(2:end);
        else
            wellseq_out =  wellseq_curr(1:end); 
            welltrigtime_out = welltrigtime_curr(1:end);
        end
        logic = zeros(size(wellseq_out)); logic(all_correct)=1;
        
        trajseq_out = 10*ones(size(wellseq_out)); %10=inbound
        trajseq_out(outbound_stidx) = 11; %11=outbound
        
        
        % Compare rewarded time with those from Output DIOS
        % Old method - only replace times where input and output triggers are close
        % -------------------------------------------------
%         output_idxs = lookup(rewtime_curr, out_trigtx); x=out_trigtx(output_idxs);
%         difft = (rewtime_curr - out_trigtx(output_idxs))./10; % in msec; /10000 for sec
%         replace = find(abs(difft)<500); % Difference is usually on the orders of tens of ms
%         rewtime_curr(replace) = x(replace);
%         % Replace in the long sequence
%         a=find(logic==1);
%         a_replace=a(replace);
%         welltrigtime_out_updated = welltrigtime_out;
%         welltrigtime_out_UPDATED(a_replace) = x(replace);
        
        
        % Replace all putative rewarded times with Output triggers. Then
        % plot with sj_plot_triggers_alltraj_2.m and check if positions are right
        % ------------------------------------------------------------------
        output_idxs = lookup(rewtime_curr, out_trigtx); % Find times in output trig close to input trig
        out_rewtimes = out_trigtx(output_idxs); % Output trigger times at reward 
        replace=find(logic==1); % Index of all the rewarded times
        %rewtime_curr = out_rewtimes; % Update reward time
        % Replace in the long sequence
        welltrigtime_out_updated = welltrigtime_out;
        welltrigtime_out_updated(replace) = out_rewtimes; % Replace rewarded times with output triggers
       
        
%         figure; hold on; 
%         plot(welltrigtime_out./10000,ones(size(welltrigtime_out)),'bo');
%         plot(out_rewtimes./10000,ones(size(out_rewtimes)),'rd');
   
        %pad = zeros(length(welltrigtime_out_updated)-length(out_rewtimes),1);
        %out_rewtimes = [out_rewtimes;pad];
        
        rewardinfo_curr = [wellseq_out, welltrigtime_out, logic, trajseq_out, welltrigtime_out_updated]; 
        rewardinfo{day}{epoch} = rewardinfo_curr;
        
    end
    
    
    
    
    eval([lowercasethree,'wellsdio = wellsdio;']); eval([lowercasethree,'rewardinfo = rewardinfo;']);
    if (directoryname(end) ~= '/')
        animdirect = [directoryname '/'];
    end
    eval(['save ',animdirect,prefix,'rewardinfo',dsz,num2str(day),' ',lowercasethree,'rewardinfo']);
    
end


