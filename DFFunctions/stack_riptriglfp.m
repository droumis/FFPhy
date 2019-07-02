function out = stack_riptriglfp(data, varargin)
% get a single matrix for the main data per animal. this makes it easier to
% work across epochs, etc
% compile results from {day}(ep).data(<ntXv>) into lfp(animal).<lfptype>(<ntXsamplesXtrial>)

% Demetris Roumis June 2019
% changed the dims on june 25 2019 to match the wavelet format.. 
% now nt x samples x rip.. will need to update the other stuff


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
    out(ian).data_dims = {'ntrode', 'sample', 'ripnum'};
    out(ian).dayeps = data(ian).epochs{1, 1};
%     if out(ian).animal == 'JZ4'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<8, :);
%     elseif out(ian).animal == 'D10'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<9, :);
%     elseif out(ian).animal == 'JZ3'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<9, :);
%     elseif out(ian).animal == 'D12'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<7, :);
%     end
    out(ian).data = {}; out(ian).numrips_perep = {};
    for ide = 1:length(out(ian).dayeps(:,1)) % day epoch
        day = out(ian).dayeps(ide,1);
        epoch = out(ian).dayeps(ide,2);
        if isempty(data(ian).output{day}(epoch).data)
            fprintf('missing data %s %d %d \n', out(ian).animal, day, epoch);
            continue
        end
        for t = 1:length(out(ian).lfptypes) % LFP type.. 
            % need to seperate out the different lfp bands so i can load
            % per type
            % ntrode x time X rip#
            tmp = data(ian).output{day}(epoch).data{t};
            if isempty(tmp)
                continue
            else
%                 out(ian).data{t}{ide} = permute(cat(3,tmp{:}), [3 2 1]);
                out(ian).data{t}{ide} = cat(3,tmp{:});
            end
            
        end
        % collect day/epoch level info
        out(ian).numrips_perep{ide,1} = length(out(ian).data{t}{ide}(1,1,:));
        out(ian).day{ide,1} = day*ones(length(out(ian).data{t}{ide}(1,1,:)), 1);
        out(ian).epoch{ide,1} = epoch*ones(length(out(ian).data{t}{ide}(1,1,:)), 1);
        out(ian).ripStartIdx{ide,1} = data(ian).output{day}(epoch).eventStartIndices;
        out(ian).ripEndIdx{ide,1} = data(ian).output{day}(epoch).eventEndIndices;
        out(ian).ripStartTime{ide,1} = data(ian).output{day}(epoch).LFPtimes(...
            data(ian).output{day}(epoch).eventStartIndices);
        out(ian).ripEndTime{ide,1} = data(ian).output{day}(epoch).LFPtimes(...
            data(ian).output{day}(epoch).eventEndIndices);
    end
    
    % collapse data across collected epochs into one matrix
    for t = 1:length(out(ian).lfptypes) % LFP type
        out(ian).data{t} = cell2mat(permute(out(ian).data{t}, [1 3 2])); % stack across time
    end
    out(ian).ntrodes = unique(cell2mat(cellfun(@(x) x(:,3), {data(ian).output{day}(epoch).index}, ...
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

