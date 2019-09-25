function out = tetinfocheck2(tetinfo,tetrodes)

% input: tetinfo struct
% prints out the area for the specified tetrodes

for d=1:length(tetinfo)
    for ep=1:length(tetinfo{d})
        if ~isempty(tetinfo{d}{ep})
            for tet=tetrodes 
                try
                   disp(sprintf('Day %d Epoch %d Tet: %d Area: %s',d,ep,tet,tetinfo{d}{ep}{tet}.area))
                catch
                   disp(sprintf('Day %d Epoch %d Tet: %d Area field not exist!',d,ep,tet))
                end
            end
        end
    end
end
                   


end

