function out = stack_riptriglfp(data, Fp, varargin)
% make single data matrix per animal
% from F(anim).output{dayep}.data{lfptype}{ripN: ntXsamp}
% into out(ian).data{lfptype}(ntrode x sample x ripple)
% exclude any events containing at least one nan in any tetrode
%
%        Comely Cactus
% 
%        _x_x__x_____x     
%       x  / x | x x  \    
%      x  x x| x |x x  x   
%      |  | |x | || |  |   
%      |  x || x x| |  x   
%     __\__x_x_|_x_x__/__  
%     \                 /  
%      `---------------'   
%       |              |    
%       \_____________/                            
%
%
%{
Notes:
- forest:bear:cactus:mushroom:leaf

FFPhy V0.1
@DR
%}
saveout = 1;
if ~isempty(varargin)
    assign(varargin)
end

out = struct;
for ian = 1:length(data) % per animal
    try
        animal = data(ian).animal{3};
    catch
        animal = data(ian).animal;
    end
    out(ian).animal = animal;
    out(ian).lfptypes = data(ian).datafilter_params.LFPtypes;
    try
        out(ian).LFPrangesHz = data(ian).datafilter_params.LFPrangesHz;
    catch
    end
    out(ian).data_dims = {'ntrode', 'sample', 'ripple'};
    dayeps = data(ian).epochs{1, 1};
    out(ian).dayeps = dayeps;
    for t = 1:length(out(ian).lfptypes) % LFP type..
        for ide = 1:size(dayeps,1) % per epoch
            day = dayeps(ide,1);
            epoch = dayeps(ide,2);
            try
                epData = data(ian).output{1,ide}.data{t};
                if isempty(epData)
                    fprintf('missing data %s %d %d\n', animal, day, epoch);
                    continue
                end
            catch
                fprintf('missing data %s %d %d\n', animal, day, epoch);
                continue
            end
            epDataTnsr = cat(3,epData{:}); % concat into [ntrode x sample x event]
            % exclude any events with a nan
            vEvents = find(squeeze(all(all(~isnan(epDataTnsr),1),2)));
            numValid = numel(vEvents);
            if max(vEvents) ~= numValid
                fprintf('D:%d E:%d LFP:%s excluding bc nan: %d of %d events\n', ...
                    day, epoch, out(ian).lfptypes{t}, max(vEvents)-length(vEvents), length(vEvents));
            end
            dataTnsr{1,ide} = epDataTnsr(:,:,vEvents);
            evStartIdx{ide,1} = data(ian).output{ide}.eventStartIndices(vEvents);
%             evEndIdx{ide,1} = data(ian).output{ide}.eventEndIndices(vEvents);
            
            % collect day/epoch level info
            numEvPerEp{ide,1} = numValid;
            evDay{ide,1} = day*ones(numValid, 1);
            evEpoch{ide,1} = epoch*ones(numValid, 1);
            evStart{ide,1} = data(ian).output{ide}.LFPtimes(evStartIdx{ide});
%             evEnd{ide,1} = data(ian).output{ide}.LFPtimes(evEndIdx{ide});
        end
%         out(ian).time{t} = data(ian).datafilter_params.time;
        out(ian).day{t} = cell2mat(evDay);
        out(ian).epoch{t} = cell2mat(evEpoch);
        out(ian).ntrodes{t} = unique(cell2mat(cellfun(@(x) x(:,3), {data(ian).output{ide}.index}, ...
            'un', 0)'), 'stable');
        
        out(ian).data{t} = cell2mat(permute(dataTnsr, [1 3 2])); % stack all rips
        out(ian).numEvPerEp{t} = cell2mat(numEvPerEp);
        out(ian).evStartIdx{t} = cell2mat(evStartIdx);
%         out(ian).evEndIdx{t} = cell2mat(evEndIdx);
        out(ian).evStart{t} = cell2mat(evStart);
%         out(ian).evEnd{t} = cell2mat(evEnd);
    end
%     % collapse data across collected epochs into one matrix
%     for t = 1:length(out(ian).lfptypes) % LFP type
%         out(ian).data{t} = cell2mat(permute(out(ian).data{t}, [1 3 2])); % stack all rips
%     end
end
if saveout
    save_data(out, Fp.paths.resultsDirectory, ['tensor_' Fp.Label]);
end

end



