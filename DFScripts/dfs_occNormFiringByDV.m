


% load in non-multiunit clusters

% log DV position in tetInfo.. idk about this.. maybe just copy from pyphy

%% 3 subplots per cell
% plot openfield firing rate maps

% plot 2d fft + canonical grid scores

% plot theta rhythmicity (how did moser's calculate intrinsic theta again?)

% plot swr response



filtfunction = 'dfa_occNormFiring';

animals = {'D12'};
days = 1;
epochfilter = '((isequal($type, ''run'')) && (isequal($environment, ''openfield'')))';
% tetfilter 
cellfilter = 'all(cellfun(''isempty'', (cellfun(@(x) strfind(x, ''mua''), $tags, ''un'', 0))))';
timefilter{1} = {'get2dstate', '(abs($velocity) >= 4)'};
iterator = 'singlecellanal';

F = createfilter('animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,...
    'excludetime', timefilter,'iterator',iterator);
F = setfilterfunction(F, filtfunction, {'spikes', 'linpos', 'pos', 'task'});
F = runfilter(F);
%%
for i = 1:size(F.output{1},2)
    try
        sp = F.output{1}(i).smoothedspikerate{1};
        image(sp)
        pause
    catch
        continue
    end
end
%%
position = [.1 .1 .8 .5];
for i = 1:size(F.output{1},2)
    try
        figure
        sp = F.output{1}(i).smoothedspikerate{1};
        image(sp)
        pause
    catch
        continue
    end
end