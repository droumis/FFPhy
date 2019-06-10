%% takes the stack of results per an,day,ep,nt,cl and returns similar but 
% MU clusters per ntrode combined

function Fout = combine_mu_clusters(F)
    animals = cellfun(@(x) x{3}, {F.animal}, 'un', 0);
    for ian = 1:numel(animals)
        Fout(ian) = F(ian); % init
        andef{ian} = animaldef(animals{ian});
        cellinfo = loaddatastruct(andef{ian}{2}, andef{ian}{3}, 'cellinfo');
        animal = animals{ian};
        fprintf('%s\n',animal);
        idata = Fout(ian).output{1};
        data_keys = cell2mat({idata.index}');
        mu_keys = evaluatefilter(cellinfo, 'isequal($tags, {''mua''})');
        mu_keys = mu_keys(find(ismember(mu_keys,data_keys, 'rows')),:);
        su_keys = data_keys(find(~ismember(data_keys, mu_keys, 'rows')),:);
%         su_keys = evaluatefilter(cellinfo, '~isequal($tags, {''mua''})');
        su_inds = find(ismember(data_keys, su_keys, 'rows'));
        % for each set of MU within the same DAY EPOCH NTRODE 
        unq_NT_keys = unique(mu_keys(:,[1:3]),'rows');
%         mu = struct;
        Fout(ian).output{1} = idata(su_inds);
        for iden = 1:length(unq_NT_keys(:,1))
            den = unq_NT_keys(iden,:);
            den_data_inds = find(ismember(data_keys(:,1:3), den, 'rows'));
            mu = idata(den_data_inds(1)); % initilize with first
            
            R = {idata(den_data_inds).psth};
            mu.psth = sum(cat(3,R{:}),3);
            
            mu.nospikes = sum([idata(den_data_inds).nospikes]);
            
            % take mean of frhist and instantFR 
            R = {idata(den_data_inds).frhist};
            mu.frhist = mean(cat(3,R{:}),3);
            
            R = {idata(den_data_inds).instantFR};
            mu.instantFR = mean(cat(3,R{:}),3);
            Fout(ian).output{1} = [Fout(ian).output{1} mu];
        end
        % add the SU and combined MU to the output
        
    end
end