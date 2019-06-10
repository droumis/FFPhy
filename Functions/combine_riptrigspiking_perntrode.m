function out = combine_riptrigspiking_perntrode(spikes, varargin)
% input data (spikes) should be struct as it is from ff

% %% takes the stack of results per an,day,ep,nt,cl and returns an,nt with 
% mu cl combined per ntrode


combine_mu = 1;
if ~isempty(varargin)
    assign(varargin)
end

% animals = cellfun(@(x) x.F.animal{3},data,'un',0');

if combine_mu
    spikes = combine_mu_clusters(spikes);
end

animals = cellfun(@(x) x{3}, {spikes.animal}, 'un', 0);
for ian = 1:numel(spikes)
    animal = animals{ian};
    andef = animaldef(animal);
    cellinfo = loaddatastruct(andef{2}, andef{3}, 'cellinfo');
    fprintf('combining per ntrode: %s\n',animal);
    idata = spikes{ian}.F.output{1};
    data_keys = cell2mat({idata.index}');
    ntrodes = unique(data_keys(:,3));
    mu_keys = evaluatefilter(cellinfo, 'isequal($tags, {''mua''})');
    mu_keys = mu_keys(find(ismember(mu_keys,data_keys, 'rows')),:);
    su_keys = data_keys(find(~ismember(data_keys, mu_keys, 'rows')),:);
    for int = 1:numel(ntrodes)
        ntrode = ntrodes(int);
        out{ian}(int).ntrode = ntrode;
        nt_mu_keys = mu_keys(mu_keys(:,3)==ntrode,:);
        nt_mu_data_inds = find(ismember(data_keys,nt_mu_keys,'rows'));
        % MU should be already combined into 1 cluster per ntrode
        out{ian}(int).mudata = idata(nt_mu_data_inds);
        % need to incorporate eplengths and psth time into the data level
        % since they are common to all ntrodes for this animal
        % this is how i was doing it previously elsewhere.. just need to
        % adapt it..
%         psth{ian}(int).eplengths{days(idy)} = cellfun(@(x) ...
%             length(x(:,1)), per_ep_mu_stack,'un',1);
%         psth{ian}(int).psth_time = idata(nt_ep_mu_inds(k)).time;
        if any(su_keys(:,3)==ntrode)
            nt_su_keys = su_keys(su_keys(:,3)==ntrode,:);
            nt_su_data_inds = find(ismember(data_keys,nt_su_keys,'rows'));
            out{ian}(int).sudata = idata(nt_su_data_inds);
        else
            out{ian}(int).sudata = struct;
        end
    end
end
end
