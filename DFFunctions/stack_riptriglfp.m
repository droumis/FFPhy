%% ---- compile results from {day}(ep).data(<ntXv>) 
% into lfp(animal).<lfptype>(<ripXtimeXnt>)

function out = stack_riptriglfp(data, varargin)
if ~isempty(varargin)
    assign(varargin)
end

out = struct;
for ian = 1:length(data)
    out(ian).animal = data{ian}.F.animal{3};
    out(ian).lfptypes = data{ian}.F.datafilter_params.LFPtypes;
    out(ian).LFPrangesHz = data{ian}.F.datafilter_params.LFPrangesHz;
    out(ian).data_dims = {'ripnum', 'time', 'ntrode'};
    out(ian).dayeps = data{ian}.F.epochs{1, 1};
%     if out(ian).animal == 'JZ4'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<8, :);
%     elseif out(ian).animal == 'JZ3'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<3, :);
%     end
    out(ian).data = {}; out(ian).numrips = []; out(ian).time = [];
    for ide = 1:length(out(ian).dayeps(:,1)) % day epoch
        day = out(ian).dayeps(ide,1);
        epoch = out(ian).dayeps(ide,2);
        for t = 1:length(out(ian).lfptypes) % LFP type
            try
                % rip# X time X ntrode
                tmp = data{ian}.F.output{day}(epoch).data{t};
                out(ian).data{t}{ide} = permute(cat(3,tmp{:}), [3 2 1]);
            catch
                continue
            end
        end
        % collect day/epoch level info
        out(ian).numrips_perep{ide,1} = length(out(ian).data{t}{ide}(:,1,1));
        out(ian).day{ide,1} = day*ones(length(out(ian).data{t}{ide}(:,1,1)), 1);
        out(ian).epoch{ide,1} = epoch*ones(length(out(ian).data{t}{ide}(:,1,1)), 1);
        out(ian).ripStartIdx{ide,1} = data{ian}.F.output{day}(epoch).eventStartIndices;
        out(ian).ripEndIdx{ide,1} = data{ian}.F.output{day}(epoch).eventEndIndices;
        out(ian).ripStartTime{ide,1} = data{ian}.F.output{day}(epoch).LFPtimes(...
            data{ian}.F.output{day}(epoch).eventStartIndices);
        out(ian).ripEndTime{ide,1} = data{ian}.F.output{day}(epoch).LFPtimes(...
            data{ian}.F.output{day}(epoch).eventEndIndices);
    end
    
    % collapse data across collected epochs into one matrix
    for t = 1:length(out(ian).lfptypes) % LFP type
        out(ian).data{t} = cell2mat(out(ian).data{t}'); % stack across time
    end
    ntrodes = unique(cell2mat(cellfun(@(x) x(:,3), {data{ian}.F.output{day}(epoch).index}, ...
        'un', 0)'), 'stable');
    out(ian).numrips_perep = cell2mat(out(ian).numrips_perep);
    out(ian).day = cell2mat(out(ian).day);
    out(ian).epoch = cell2mat(out(ian).epoch);
    out(ian).ripStartIdx = cell2mat(out(ian).ripStartIdx);
    out(ian).ripEndIdx = cell2mat(out(ian).ripEndIdx);
    out(ian).ripStartTime = cell2mat(out(ian).ripStartTime);
    out(ian).ripEndTime = cell2mat(out(ian).ripEndTime);
end

