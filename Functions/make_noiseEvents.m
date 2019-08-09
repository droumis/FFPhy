

function out = make_noiseEvents(noizevents, lfpstack, varargin)
%{
Inputs
noizevents = mat of indices into lfpstack:.day.epoch.ripStartTimes
%}
saveout = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
animals = {lfpstack.animal};
for ian = 1:length(animals)
    animal = animals{ian};
    lfpanidx = find(strcmp({lfpstack.animal}, animal));
    nzanidx = find(strcmp({noizevents.animal}, animal));
    out(ian).animal = animal;
    out(ian).dims = {'event', {'day', 'epoch', 'time'}};
    r = noizevents(nzanidx).ripnums;
    out(ian).events = [lfpstack(lfpanidx).day(r) lfpstack(lfpanidx).epoch(r) ...
        lfpstack(lfpanidx).ripStartTime(r)];
    days = unique(out(ian).events(:,1))';
    for d = days
        dr = out(ian).events(:,1) == d;
        eps = unique(out(ian).events(dr,2))';
        for e = eps
            deridx = ismember(out(ian).events(:,1:2), [d e], 'rows');
            if out(ian).events(deridx,2) ~= e
                error('bad ep');
            elseif out(ian).events(deridx,1) ~= d
                error('bad day')
            end
            out(ian).perepoch{d}{e}{1}.starttime = out(ian).events(deridx,3);
        end
    end
    if saveout
        andef = animaldef(animal);
        save_data(out(ian), andef{2}, 'noiseEvents');
    end
end