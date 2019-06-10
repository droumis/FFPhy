
function [wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt_lin(directoryname,prefix,days,epochs, dionos, varargin)
% [wellsdio] =
% sj_findwellsfromdio1('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',[1:3],[2 4],[26 27 28]);
% [wellsdio] = sj_findwellsfromdio1('/data25/sjadhav/RippleInterruption/REd_direct','REd',[1:2],[2 4],[22 23 24]);
% [wellsdio] = sj_findwellsfromdio1('/data25/sjadhav/HPExpt/HPb_direct','HPb',[1],[4 6],[24 23 22]);
% sj_findwellsfromdio1('/data25/sjadhav/HPExpt/HPa_direct','HPa',[2],[2 4],[24 23 22]);
%% From DIO, gets well start and end for all unique trajectories
format long
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
    
    %% EGYPT LINEAR (need to verify)
    %usedios = [24,23];
    %outdios=[8,7];
    
    %% DAVE LINEAR
    usedios = [23,24];      % INPUTS: well 1, well2
    outdios=[7,8];        % OUTPUTS: well 1, well 2
    
    
    
    for epoch=epochs,
        
        %Init
        well_start=[]; well_end=[]; well_startend=[]; wellseq_curr=[]; welltrigtime_curr=[]; tmpwells=[];
        
        DIO1=DIO{day}{epoch}{usedios(1)}; nones = length(DIO1.pulsetimes); % Center well 1: Bits 8/24 in
        DIO2=DIO{day}{epoch}{usedios(2)}; nzeros = length(DIO2.pulsetimes); % Well 0: Bits 9/25 in
%         DIO3=DIO{day}{epoch}{usedios(3)}; ntwos = length(DIO3.pulsetimes); % Well 2: Bits 10/26 in
%         
        % Output wells
        DIOout1=DIO{day}{epoch}{outdios(1)}; nones_out = length(DIOout1.pulsetimes); % Center well 1: Bits 8 out/24
        DIOout2=DIO{day}{epoch}{outdios(2)}; nzeros_out = length(DIOout2.pulsetimes); % Well 0: Bits 9 out/25
%         DIOout3=DIO{day}{epoch}{outdios(3)}; ntwos_out = length(DIOout3.pulsetimes); % Well 2: Bits 10 out/26
        out_trigt = [DIOout1.pulsetimes(:,1);DIOout2.pulsetimes(:,1)];
        [out_trigtx, out_ord] = sort(out_trigt);
        out_wells = [1*ones(nones_out,1);0*ones(nzeros_out,1)];
        out_rewwells = out_wells(out_ord);
        
        wells = [1*ones(nones,1);0*ones(nzeros,1)];
        trigt = [DIO1.pulsetimes(:,1);DIO2.pulsetimes(:,1)];
        [tfx,ti] = sort(trigt);
        wfx=wells(ti);
        swit=find(diff(wfx)~=0)+1;
        well_start=[wfx(1);wfx(swit(1:end-1))];
        well_end=[wfx(swit)];
        well_startend=[well_start,well_end];
        % For saving
        wellsdio{day}{epoch}=well_startend;
        
        % NEW
        % ----
        wellseq_curr=[wfx(1);wfx(swit)];
        welltrigtime_curr=[tfx(1);tfx(swit)];
        
        % IF First Well is 1, then first outbound was captured. Keep it aside temporarily
         % Remove first outbound if it exists. Put back later
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
            
        %% Find outbound (rat starts at well 8 (1) on linear track)
        outbound_stidx = find(tmpwells(:,1)==1);
        outbound_wellstend = tmpwells(outbound_stidx,:);
        outbound_logic = zeros(length(outbound_wellstend),1);
        corridx = find( (tmpwells(outbound_stidx,2)~=tmpwells(outbound_stidx,1))...
            & (tmpwells(outbound_stidx,2)~=1) );                                    %Finds trips where end well does not equal start well
        correct_out = outbound_stidx(corridx);
        wrong_out = setdiff(outbound_stidx,correct_out);
        %% Find inbound (rat starts at well 9 (0) on linear track)
        inbound_stidx = find(tmpwells(:,1)~=1);
        inbound_wellstend = tmpwells(inbound_stidx,:);
        % Find out which inbound is correct, and return the logic
        inbound_logic = zeros(length(inbound_wellstend),1);
        corridx = find( (tmpwells(inbound_stidx,2)~=tmpwells(inbound_stidx,1))...
            & (tmpwells(inbound_stidx,2)==1) );
        correct_in = inbound_stidx(corridx);
        wrong_in = setdiff(inbound_stidx,correct_in);
        %inbound_logic(corr) = 1;
        
        % IMP - Put first Outbound Back In If It Was Removed - Have to Push all indexes By 1
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
                
        all_correct = sort([correct_out; correct_in]);
        all_wrong = sort([wrong_out; wrong_in]);
        allwell_curr = sort([all_correct;all_wrong]);
        rewtime_curr = welltrigtime_curr(all_correct+1);
        norewtime_curr = welltrigtime_curr(all_wrong+1);
        rewardedwell_seq = well_startend(all_correct,2); unrewardedwell_seq = well_startend(all_wrong,2);
        wellseq_reconstructed = wellseq_curr(allwell_curr+1);
        
        % Now prepare outputs
        wellseq_out =  wellseq_curr(2:end); 
        welltrigtime_out = welltrigtime_curr(2:end);
        logic = zeros(size(wellseq_out)); logic(all_correct)=1;
        
        trajseq_out = 10*ones(size(wellseq_out)); %10=inbound
        trajseq_out(outbound_stidx) = 11;
        
        
        % Compare rewarded time with those from Output DIOS
        output_idxs = lookup(rewtime_curr, out_trigtx); x=out_trigtx(output_idxs);
        difft = (rewtime_curr - out_trigtx(output_idxs))./10; % in msec; /10000 for sec
        replace = find(abs(difft)<600); % Difference is usually on the orders of tens of ms
        rewtime_curr(replace) = x(replace);
        % Replace in the long sequence
        a=find(logic==1);
        a_replace=a(replace);
        welltrigtime_out_updated = welltrigtime_out;
        welltrigtime_out(a_replace) = x(replace);
        
        
        rewardinfo_curr = [wellseq_out, welltrigtime_out, logic, trajseq_out, welltrigtime_out_updated]; 
        rewardinfo{day}{epoch} = rewardinfo_curr;
        
    end
    %wellseq_out = ending well of trajectory (column 2 of wellsdio)
    %welltrigtime_out = input trigger time at that well
    %logic = 1 is correct, 0 is incorrect
    %trajsec_out = 10 is inbound trial, 11 is outbound trial; does not
    %    reflect whether trial is correct, only location of well
    %welltrigtime_out_updated = output trigger time at that well (only
    %    updated if correct trial); NOTE: for trials with 500 ms delay between input and output, make sure abs(difft)>500, to 
    %    include output times that are within something greater than 500 ms from the input time 
    
    
    
    
    eval([lowercasethree,'wellsdio = wellsdio;']); eval([lowercasethree,'rewardinfo = rewardinfo;']);
    if (directoryname(end) ~= '/')
        animdirect = [directoryname '/'];
    end
    eval(['save ',animdirect,prefix,'rewardinfo',dsz,num2str(day),' ',lowercasethree,'rewardinfo']);
    
end


