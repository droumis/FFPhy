function out = tetinfocheck(tetinfo)

% find the Reference tetrode

% input: tetinfo struct
% prints out the reference tetrode for every epoch

for d=1:length(tetinfo)
    for ep=1:length(tetinfo{d})
        if ~isempty(tetinfo{d}{ep})
            for tet=1:length(tetinfo{d}{ep}) 
                try
                   area = tetinfo{d}{ep}{tet}.area;
                   disp(sprintf('Day %d Epoch %d Tet: %d Area: %s',d,ep,tet,area))
                   %if (area == 'Reference')
                   %    disp(sprintf('Day %d Epoch %d Ref: %d',d,ep,tet))
                   %end
                catch
                end
            end
        end
    end
end
                   


end

