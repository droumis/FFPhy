function out = stack_riptriglfp(data, Fp, varargin)
% make single data matrix per animal

% compile results from F(anim).output{dayep}.data{lfptype}{ripN: ntXsamp}
% into out(ian).data{lfptype}(ntrode x sample x ripple)

% Demetris Roumis June 2019
saveout = 1;
if ~isempty(varargin)
    assign(varargin)
end

out = struct;
for ian = 1:length(data)
    out(ian).animal = data(ian).animal{3};
    out(ian).lfptypes = data(ian).datafilter_params.LFPtypes;
    try
        out(ian).LFPrangesHz = data(ian).datafilter_params.LFPrangesHz;
    catch
    end    
    out(ian).data_dims = {'ntrode', 'sample', 'ripple'};
    out(ian).dayeps = data(ian).epochs{1, 1};
    out(ian).data = {};
    out(ian).numrips_perep = {};
    for ide = 1:length(out(ian).dayeps(:,1)) % per epoch
        day = out(ian).dayeps(ide,1);
        epoch = out(ian).dayeps(ide,2);
        try
        if isempty(data(ian).output{ide}.data)
            fprintf('missing data %s %d %d\n', out(ian).animal, day, epoch);
            continue
        end
        catch
            fprintf('missing data %s %d %d\n', out(ian).animal, day, epoch);
            continue
        end
        for t = 1:length(out(ian).lfptypes) % LFP type.. 
            % TODO: need to seperate out the different lfp bands so i can load per type
            % ntrode x sample X ripple
            tmp = data(ian).output{1,ide}.data{t};
            if isempty(tmp)
                continue
            else
                out(ian).data{t}{1,ide} = cat(3,tmp{:}); % concat ripples into 3rd dim
            end
        end
        % collect day/epoch level info
        out(ian).numrips_perep{ide,1} = length(out(ian).data{t}{ide}(1,1,:));
        out(ian).day{ide,1} = day*ones(length(out(ian).data{t}{ide}(1,1,:)), 1);
        out(ian).epoch{ide,1} = epoch*ones(length(out(ian).data{t}{ide}(1,1,:)), 1);
        out(ian).ripStartIdx{ide,1} = data(ian).output{ide}.eventStartIndices;
        out(ian).ripEndIdx{ide,1} = data(ian).output{ide}.eventEndIndices;
        out(ian).ripStartTime{ide,1} = data(ian).output{ide}.LFPtimes(...
            data(ian).output{ide}.eventStartIndices);
        out(ian).ripEndTime{ide,1} = data(ian).output{ide}.LFPtimes(...
            data(ian).output{ide}.eventEndIndices);
    end
    % collapse data across collected epochs into one matrix
    for t = 1:length(out(ian).lfptypes) % LFP type
        out(ian).data{t} = cell2mat(permute(out(ian).data{t}, [1 3 2])); % stack all rips
    end
    out(ian).ntrodes = unique(cell2mat(cellfun(@(x) x(:,3), {data(ian).output{ide}.index}, ...
        'un', 0)'), 'stable');
    out(ian).time = data(ian).datafilter_params.time;
    out(ian).numrips_perep = cell2mat(out(ian).numrips_perep);
    out(ian).day = cell2mat(out(ian).day);
    out(ian).epoch = cell2mat(out(ian).epoch);
    out(ian).ripStartIdx = cell2mat(out(ian).ripStartIdx);
    out(ian).ripEndIdx = cell2mat(out(ian).ripEndIdx);
    out(ian).ripStartTime = cell2mat(out(ian).ripStartTime);
    out(ian).ripEndTime = cell2mat(out(ian).ripEndTime);
end
if saveout
    save_data(out, Fp.paths.resultsDirectory,['riptriglfpstack_',Fp.epochEnvironment]);
end

end


