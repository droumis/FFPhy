%% takes the stack of results per an,day,ep,nt,cl and returns an,nt with cl
% combined and epochs vertstacked

function out = combine_riptrigspiking_perntrode(data, varargin)

combine_mu = 1;
if ~isempty(varargin)
    assign(varargin)
end
animals = cellfun(@(x) x.F.animal{3},data,'un',0');
for ian = 1:numel(animals)
    if combine_mu
        data{ian}.F = combine_mu_clusters(data{ian}.F);
    end
    andef = animaldef(animals{ian});
    cellinfo = loaddatastruct(andef{2}, andef{3}, 'cellinfo');
    animal = animals{ian};
    fprintf('combining per ntrode: %s\n',animal);
    idata = data{ian}.F.output{1};
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

%             % don't ep stack anymore
%
%             nt_uniq_dayeps = unique(nt_data_keys(:,[1 2]),'rows');
% %             days = unique(nt_data_keys(:,1));
%
%             psth{ian}(int).days = days;
%             for ide = 1:length(nt_uniq_dayeps(:,1)) % day
%                 day = nt_uniq_dayeps(ide,:)
%                 per_ep_mu_stack = {};
%                 nt_day_keys = nt_data_keys(nt_data_keys(:,1)==days(ide,:),:);
%
%                 % find and combine cluster results and across eps
%                 nt_mu_keys = nt_day_keys(find(ismember(nt_day_keys, ...
%                     mu_keys, 'rows')),:);
%                 nt_su_keys = nt_day_keys(find(ismember(nt_day_keys, ...
%                     su_keys, 'rows')),:);
%                 if ~isempty(nt_mu_keys)
%                     epochs = unique(nt_mu_keys(:,2));
%                     for iep = 1:numel(epochs) % epoch
%                         per_ep_mu_stack{iep,1} = 0;
%                         epoch = epochs(iep);
%                         nt_ep_mu_keys = nt_mu_keys(nt_mu_keys(:,2)==epoch,:);
%                         nt_ep_mu_inds = find(ismember(data_keys, nt_ep_mu_keys, ...
%                             'rows'));
%                         for k = 1:numel(nt_ep_mu_inds) % MU cluster
%                             try
%                                 per_ep_mu_stack{iep,1} = per_ep_mu_stack{iep,1} + ...
%                                     idata(nt_ep_mu_inds(k)).psth;
%                             catch
%                                 fprintf('could not add mu spikes from %d %d %d %d, skipping\n',nt_ep_mu_keys(k,:));
%                                 continue
%                             end
%                         end
%                     end
%                     psth{ian}(int).eplengths{days(ide)} = cellfun(@(x) ...
%                         length(x(:,1)), per_ep_mu_stack,'un',1);
%                     psth{ian}(int).psth_time = idata(nt_ep_mu_inds(k)).time;
%                 end
%                 psth{ian}(int).mucluster{days(ide)} = {cell2mat(per_ep_mu_stack)};
%                 if ~isempty(nt_su_keys)
%                     sus = unique(nt_su_keys(:,4));
%                     for isu = 1:numel(sus) % SU clusters
%                         sui_keys = nt_su_keys(ismember(nt_su_keys(:,4),sus(isu),'rows'),:);
%                         nt_ep_su_inds = find(ismember(data_keys, sui_keys, ...
%                             'rows'));
%                         psth{ian}(int).suclusters{days(ide)}{isu} = cell2mat({idata(nt_ep_su_inds).psth}');
%                         psth{ian}(int).sucluster_id{days(ide)}{isu} = sus(isu);
%                     end
%                 end
%             end
%         end
%     end
% end