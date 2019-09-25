
function [wellsdio, rewardinfo] = sj_findwellsfromdio1_lintrack(directoryname,prefix,days,epochs, dionos, varargin)

% From sj_findwellsfromdio1. Adapt the w-track to lintrack - which is easy. Only look at outputs
% sj_findwellsfromdio1_lintrack('/data25/sjadhav/HPExpt/HPb_direct','HPb',[1],[2],[6 7]);

% For Wtrack, it was:
% [wellsdio] = sj_findwellsfromdio1('/data25/sjadhav/HPExpt/HPb_direct','HPb',[1],[4 6],[24 23 22]);

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
 
    outdios=dionos;
    
    for epoch=epochs,
        
        %Init
        well_startend=[]; wellseq_curr=[]; welltrigtime_curr=[]; tmpwells=[];
          
        % Output wells
        DIOout1=DIO{day}{epoch}{outdios(1)}; nones_out = length(DIOout1.pulsetimes); %  well 1
        DIOout2=DIO{day}{epoch}{outdios(2)}; nzeros_out = length(DIOout2.pulsetimes); % Well 0:
        out_trigt = [DIOout1.pulsetimes(:,1);DIOout2.pulsetimes(:,1)];
        [out_trigtx, out_ord] = sort(out_trigt);
        out_wells = [1*ones(nones_out,1);0*ones(nzeros_out,1)];
        out_rewwells = out_wells(out_ord);
        logic = ones(size(out_wells));
        well_startend = [out_rewwells(1:end-1),out_rewwells(2:end)];
        % For saving
        wellsdio{day}{epoch}=well_startend;
        
        welltrigtime_out = out_trigtx; welltrigtime_out_updated = welltrigtime_out;
        trajseq_out = out_rewwells;
        rewardinfo_curr = [out_rewwells, welltrigtime_out, logic, trajseq_out, welltrigtime_out_updated]; 
        rewardinfo{day}{epoch} = rewardinfo_curr;
      
    end
    
    
    
    
    eval([lowercasethree,'wellsdio = wellsdio;']); eval([lowercasethree,'rewardinfo = rewardinfo;']);
    if (directoryname(end) ~= '/')
        animdirect = [directoryname '/'];
    end
    eval(['save ',animdirect,prefix,'rewardinfo_lintrack',dsz,num2str(day),' ',lowercasethree,'rewardinfo']);
    
end


