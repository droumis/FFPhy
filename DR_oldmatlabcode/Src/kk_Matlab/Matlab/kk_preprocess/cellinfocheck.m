function out = cellinfocheck(cellinfo)

% find the Reference tetrode

% input: tetinfo struct
% prints out the reference tetrode for every epoch

for d=1:length(cellinfo)
    for ep=1:length(cellinfo{d})
        if ~isempty(cellinfo{d}{ep})
            for tet=1:length(cellinfo{d}{ep})
                for cell=1:length(cellinfo{d}{ep}{tet})
                    try
                        area = cellinfo{d}{ep}{tet}{cell}.area;
                        type = cellinfo{d}{ep}{tet}{cell}.type;
                        disp(sprintf('Day %d Epoch %d Tet: %d Area: %s Type: %s',d,ep,tet,area,type))
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



end

