function Fout = stack_riptrigspiking(F, varargin)
% Z stack mu and su matrices for each dayep

% Demetris Roumis June 2019

if ~isempty(varargin)
    assign(varargin)
end

% function Fout = stack_riptrigspiking(F)
animals = cellfun(@(x) x{3}, {F.animal}, 'un', 0);
F = combine_mu_clusters(F);
for ian = 1:numel(animals)
    %         Fout(ian) = F(ian); % init
    andef{ian} = animaldef(animals{ian});
    cellinfo = loaddatastruct(andef{ian}{2}, andef{ian}{3}, 'cellinfo');
    animal = animals{ian};
    fprintf('stacking %s\n',animal);
    idata = F(ian).output{1};
    data_keys = cell2mat({idata.index}');

    mu_keys = evaluatefilter(cellinfo, 'isequal($tags, {''mua''})');
    mu_keys = mu_keys(ismember(mu_keys,data_keys, 'rows'),:);
    su_keys = data_keys(~ismember(data_keys, mu_keys, 'rows'),:);
%     su_keys = evaluatefilter(cellinfo, '~isequal($tags, {''mua''})');
    
    ntrodes = unique(data_keys(:,3));
    Fout(ian).ntrodes = ntrodes;
    unq_DE_keys = unique(data_keys(:,[1:2]),'rows');
    if animal == 'JZ1'
        unq_DE_keys = unq_DE_keys(unq_DE_keys(:,1)<7,:);
    elseif animal == 'D10'
        unq_DE_keys = unq_DE_keys(unq_DE_keys(:,1)<9,:);
    elseif animal == 'JZ3'
        unq_DE_keys = unq_DE_keys(unq_DE_keys(:,1)<9,:);
    elseif animal == 'D12'
        unq_DE_keys = unq_DE_keys(unq_DE_keys(:,1)<7,:);
    elseif animal == 'JZ4'
        unq_DE_keys = unq_DE_keys(unq_DE_keys(:,1)<7,:);
    end
    for ide = 1:length(unq_DE_keys(:,1))
        de = unq_DE_keys(ide,:);
        
        su_de_keys = su_keys(ismember(su_keys(:,1:2), de, 'rows'),:);
        su_de_inds = find(ismember(data_keys, su_de_keys, 'rows'));
        mu_de_keys = mu_keys(find(ismember(mu_keys(:,1:2), de, 'rows')),:);
        mu_de_inds = find(ismember(data_keys,mu_de_keys, 'rows'));
        Fout(ian).eventtimes{de(1)}{de(2)} = idata(mu_de_inds(1)).eventtags(:,2);
        if ~isempty(mu_de_inds) % combine the multi unit clusters per ntrode
            R = {idata(mu_de_inds).psth};
            Fout(ian).muspikes{de(1)}{de(2)} = cat(3,R{:});
            R = {idata(mu_de_inds).frhist};
            Fout(ian).mufrhist{de(1)}{de(2)} = cat(3,R{:});
            R = {idata(mu_de_inds).instantFR};
            Fout(ian).muinstantFR{de(1)}{de(2)} = cat(3,R{:});
            Fout(ian).mu_de_keys{de(1)}{de(2)} = mu_de_keys;
        end
        if ~isempty(su_de_inds)
            R = {idata(su_de_inds).psth};
            Fout(ian).suspikes{de(1)}{de(2)} = cat(3,R{:});
            Fout(ian).sukeys{de(1)}{de(2)} = su_de_keys;
        end
    end
end
end


%     %% multiunit
%     for nti = 1:length(ntrodes)
%         ntrode = ntrodes(nti);
% %         Fout(ian).data{ntrode} = [];
%         nt_data_keys = data_keys(find(data_keys(:,3)==ntrode),:);
%         nt_mu_inds = find(ismember(nt_data_keys(:,1:3), mu_keys(:,1:3), 'rows'));
%         nt_mu_keys = data_keys(nt_mu_inds,:);
%
%         unq_NT_keys = unique(nt_data_keys(:,[1:3]),'rows');
%         for iden = 1:length(unq_NT_keys(:,1))
%             den = unq_NT_keys(iden,1:3);
%             den_mu_inds = find(ismember(nt_mu_keys(:,1:3), den, 'rows'));
%             % becase mu has the same dimensionality as lfp, save as single
%             % mat, just like the lfp.
%             if ~isempty(den_mu_inds)
%                 R = {idata(den_mu_inds).psth};
%                 Fout(ian).muspikes{ntrode}{iden,1} = nansum(cat(3,R{:}),3);
%                 R = {idata(den_mu_inds).frhist};
%                 Fout(ian).mufrhist{ntrode}{iden,1} = nanmean(cat(3,R{:}),3);
%                 R = {idata(den_mu_inds).instantFR};
%                 Fout(ian).muinstantFR{ntrode}{iden,1} = nanmean(cat(3,R{:}),3);
%                 de_numrips = length(R{1}(:,1));
%                 Fout(ian).muday{iden,1} = repmat(den(1), de_numrips, 1);
%                 Fout(ian).muepoch{iden,1} = repmat(den(2), de_numrips, 1);
%             end
% %             sp = idata(den_data_inds(1)); % initilize with first
%         end
%         % stack across dayeps for this ntrode
%         try
%             Fout(ian).muspikes{ntrode} = cell2mat(Fout(ian).muspikes{ntrode});
%             Fout(ian).mufrhist{ntrode} = cell2mat(Fout(ian).mufrhist{ntrode});
%             Fout(ian).muinstantFR{ntrode} = cell2mat(Fout(ian).muinstantFR{ntrode});
%         catch
%             fprintf('no mu for %s nt%d\n', animal, ntrode);
%             continue
%         end
%     end
%     % stack Z ntrodes
%     Fout(ian).muspikes = cell2mat(permute(Fout(ian).muspikes, [3 1 2]));
%     Fout(ian).mufrhist = cell2mat(permute(Fout(ian).mufrhist, [3 1 2]));
%     Fout(ian).muinstantFR = cell2mat(permute(Fout(ian).muinstantFR, [3 1 2]));
%     Fout(ian).muday = cell2mat(Fout(ian).muday);
%     Fout(ian).muepoch = cell2mat(Fout(ian).muepoch);
%
%% single unit. collect per dayep since there are variable # clusters

%
% % spikes = combine_mu_clusters(spikes);
%
% out = struct;
% for ian = 1:length(spikes)
%     animal = spikes(ian).animal{3};
%
%     fprintf('stacking: %s\n',animal);
%     iandata = spikes(ian).output{1};
%     data_keys = cell2mat({iandata.index}');
%     ntrodes = unique(data_keys(:,3));
%     out(ian).ntrodes = ntrodes;
%     out(ian).dayeps = unique(data_keys(:,[1 2]), 'rows');
%     out(ian).animal = data(ian).animal{3};
%     out(ian).data_dims = {'ripnum', 'time', 'ntrode'};
%
%     for int = 1:numel(ntrodes)
%         ntrode = ntrodes(int);
%         for ide = 1:length(out(ian).dayeps(:,1)) % day epoch
%             day = out(ian).dayeps(ide,1);
%             epoch = out(ian).dayeps(ide,2);
%             tmp = iandata.data{t};
%             out(ian).data{t}{ide} = permute(cat(3,tmp{:}), [3 2 1]);
%                 store{ntrode}{dayep} =
%             end
%         store{ntrode} = cell2mat(store{ntrode});
%
%         %hacking for now
%     if out(ian).animal == 'JZ4'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<8, :);
%     elseif out(ian).animal == 'D10'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<9, :);
%     elseif out(ian).animal == 'D12'
%         out(ian).dayeps = out(ian).dayeps(out(ian).dayeps(:,1)<9, :);
%     end
%
%     out(ian).data = {}; out(ian).numrips_perep = {};
%
%
% %         if isempty(data(ian).output{day}(epoch).data)
% %             fprintf('missing data %s %d %d \n', out(ian).animal, day, epoch);
% %             continue
% %         end
%         for t = 1:length(out(ian).lfptypes) % LFP type
%             % rip# X time X ntrode
%                 tmp = data(ian).output{day}(epoch).data{t};
%             out(ian).data{t}{ide} = permute(cat(3,tmp{:}), [3 2 1]);
%         end
%         % collect day/epoch level info
%         out(ian).numrips_perep{ide,1} = length(out(ian).data{t}{ide}(:,1,1));
%         out(ian).day{ide,1} = day*ones(length(out(ian).data{t}{ide}(:,1,1)), 1);
%         out(ian).epoch{ide,1} = epoch*ones(length(out(ian).data{t}{ide}(:,1,1)), 1);
%         out(ian).ripStartIdx{ide,1} = data(ian).output{day}(epoch).eventStartIndices;
%         out(ian).ripEndIdx{ide,1} = data(ian).output{day}(epoch).eventEndIndices;
%         out(ian).ripStartTime{ide,1} = data(ian).output{day}(epoch).LFPtimes(...
%             data(ian).output{day}(epoch).eventStartIndices);
%         out(ian).ripEndTime{ide,1} = data(ian).output{day}(epoch).LFPtimes(...
%             data(ian).output{day}(epoch).eventEndIndices);
%     end
%
%     % collapse data across collected epochs into one matrix
%     for t = 1:length(out(ian).lfptypes) % LFP type
%         out(ian).data{t} = cell2mat(out(ian).data{t}'); % stack across time
%     end
%     out(ian).ntrodes = unique(cell2mat(cellfun(@(x) x(:,3), {data(ian).output{day}(epoch).index}, ...
%         'un', 0)'), 'stable');
%     out(ian).time = data(ian).datafilter_params.time;
%     out(ian).numrips_perep = cell2mat(out(ian).numrips_perep);
%     out(ian).day = cell2mat(out(ian).day);
%     out(ian).epoch = cell2mat(out(ian).epoch);
%     out(ian).ripStartIdx = cell2mat(out(ian).ripStartIdx);
%     out(ian).ripEndIdx = cell2mat(out(ian).ripEndIdx);
%     out(ian).ripStartTime = cell2mat(out(ian).ripStartTime);
%     out(ian).ripEndTime = cell2mat(out(ian).ripEndTime);
% end
%
