%% plots selected tetrodes' EEG side-by-side for, random start times

function out = plotEEG(eegstruct,day,epoch,tetrodes,no_windows,duration)

%% tetrodes is a vector of tetrodes
%% duration is in seconds

num_tetrodes=length(tetrodes);
epoch_starttime=eegstruct{day}{epoch}{1}.starttime;
no_samp=length(eegstruct{day}{epoch}{1}.data);
samprate=eegstruct{day}{epoch}{1}.samprate;
duration=floor(duration*samprate);

for w=1:no_windows
    figure
    start=floor(rand*no_samp);
    hold on
    for t=tetrodes
        if t==10 || t==12        %% CA3
%             plot(1:(duration+1), ... 
%                 eegstruct{day}{epoch}{t}.data(start:(start+duration)), ...
%                 'r','LineWidth',1);
        elseif t==8 || t==9        %% CA2
            plot(1:(duration+1), ... 
                eegstruct{day}{epoch}{t}.data(start:(start+duration)), ...
                'g','LineWidth',1);
        else
            plot(1:(duration+1), ... 
                eegstruct{day}{epoch}{t}.data(start:(start+duration)), ...
                'k','LineWidth',1);
        end
    end
    
end


end



