function output = calcsinglecellreactivationstats(g,decodefilter)

nov = [];
for an = 1:length(g)
    tmpout = calcfiringratestats2(g(an),decodefilter(an));
    nov = [nov tmpout.novel];
end


output.index = []; output.typen = []; output.types = [];
output.peakf = []; output.peakn = []; output.peaks = [];
output.fmean = []; output.nmean = []; output.smean = [];
output.activationf = []; output.activationn = []; output.activations = [];
output.swrratef = []; output.swrraten = []; output.swrrates = [];
output.locf = []; output.locn = []; output.locs = [];

for i = 1:length(nov)
    if ~isempty(nov(i).peaks)
        %Calculate the correlation between novel and familiar
        output.index = [output.index; nov(i).index];
        
        output.typen = [output.typen; nov(i).typen];
        output.types = [output.types; nov(i).types];
        
        output.peakf = [output.peakf; nov(i).peakf];
        output.peakn = [output.peakn; nov(i).peakn];
        output.peaks = [output.peaks; nov(i).peaks];
        
        output.fmean = [output.fmean; nov(i).fmean];
        output.nmean = [output.nmean; nov(i).nmean];
        output.smean = [output.smean; nov(i).smean];
        
        output.activationf = [output.activationf; nov(i).activationf];
        output.activationn = [output.activationn; nov(i).activationn];
        output.activations = [output.activations; nov(i).activations];
        
        output.swrratef = [output.swrratef; nov(i).swrratef];
        output.swrraten = [output.swrraten; nov(i).swrraten];
        output.swrrates = [output.swrrates; nov(i).swrrates];
        
        output.locf = [output.locf; nov(i).swrlocf];
        output.locn = [output.locn; nov(i).swrlocn];
        output.locs = [output.locs; nov(i).swrlocs];
       
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = calcfiringratestats2(f,decodefilter)

%Initialize output
out = [];

% Determine the index for familiar, novel, and super novel sessions
fam = [];   nov = [];   sup = [];
for i = 1:length(f.output{1});	tmp = f.output{1}(i).index; tmp(:,5) = i; fam = [fam; tmp];	end

for i = 1:length(f.output{2});
    if ~isempty(f.output{2}(i).familiarsession)
        tmp = f.output{2}(i).index; tmp(:,5) = i;   
        tmp(:,6) = f.output{2}(i).familiarsession; nov = [nov; tmp];
    end
end

if length(f.output) == 3
    for i = 1:length(f.output{3});
        if ~isempty(f.output{3}(i).familiarsession)
            tmp = f.output{3}(i).index; tmp(:,5) = i;
            tmp(:,6) = f.output{3}(i).familiarsession;
            tmp(:,7) = f.output{3}(i).novelsession; sup = [sup; tmp];
        end
    end
end
clear tmp


if ~isempty(fam)
    % Load task structure
    eval(['task',' = loaddatastruct(f.animal{2}, f.animal{3}, ''task'', unique(fam(:,1)));']);

    % For each super novel cell, calculate the peak and average difference in firing rate
    for cellindex = 1:size(sup,1)
        
        supcell = sup(cellindex,:);
        ind = rowfind(supcell(1:4),f.output{3}(supcell(5)).index);
        if ~isnan(f.output{3}(supcell(5)).quadrants(ind))
            tmp.index = supcell([1 3 4]);
            
            %Determine the quadrant and quadrant type
            tmp.types = f.output{3}(supcell(5)).type(ind);
            
            %Peak Rate
            tmp.peaks = f.output{3}(supcell(5)).peak(ind);

            %Mean Rate
            tmp.smean = f.output{3}(supcell(5)).meanrate(ind);

            %Activation probability and SWR rate
            if ~isempty(decodefilter.output{3}(supcell(5)).activationprobability)
                tmp.activations = decodefilter.output{3}(supcell(5)).activationprobability(ind);
                tmp.swrrates = decodefilter.output{3}(supcell(5)).SWRrate(ind);
                tmp.swrlocs = decodefilter.output{3}(supcell(5)).SWRlocation(ind,[2 3]);
            else
                tmp.activations = NaN;
                tmp.swrrates = NaN;
                tmp.swrlocs = nan(1,2);
            end

            %Ignore if the rate is not defined for super novel cell
            if ~isempty(tmp.peaks)
                %Find the matching novel cell
                novcell = nov(find(rowfind(nov(:,1:4),supcell([1 7 3 4]))),:);

                if ~isempty(novcell)
                    indn =rowfind(novcell(1:4),f.output{2}(novcell(5)).index);
                    if ~all(isnan(f.output{2}(novcell(5)).peak(indn)));

                        %Type
                        tmp.typen = f.output{2}(novcell(5)).type(indn);
                        %Peak rate
                        tmp.peakn = f.output{2}(novcell(5)).peak(indn);

                        %Mean rate
                        tmp.nmean = f.output{2}(novcell(5)).meanrate(indn);

                        %Activation probability and SWR rate
                        if ~isempty(decodefilter.output{2}(novcell(5)).activationprobability)
                            tmp.activationn = decodefilter.output{2}(novcell(5)).activationprobability(indn);
                            tmp.swrraten = decodefilter.output{2}(novcell(5)).SWRrate(indn);
                            tmp.swrlocn = decodefilter.output{2}(novcell(5)).SWRlocation(indn,[2 3]);
                        else
                            tmp.activationn = NaN;
                            tmp.swrraten = NaN;
                            tmp.swrlocn = nan(1,2);
                        end
                    else
                        tmp.typen = nan(1,1); tmp.peakn = nan(1,1); tmp.nmean = nan(1,1); tmp.activationn = nan(1,1); tmp.swrraten = nan(1,1); tmp.swrlocn = nan(1,2);
                    end
                else
                    tmp.typen = nan(1,1); tmp.peakn = nan(1,1); tmp.nmean = nan(1,1); tmp.activationn = nan(1,1); tmp.swrraten = nan(1,1); tmp.swrlocn = nan(1,2);
                end

                %Find the matching familiar cell 
                famcell = fam(find(rowfind(fam(:,1:4),supcell([1 6 3 4]))),:);
                if ~isempty(famcell)
                    indf = rowfind(famcell(1:4),f.output{1}(famcell(5)).index);
                    if ~isnan(f.output{1}(famcell(5)).peak(indf));
                        %Peak Rate
                        tmp.peakf = f.output{1}(famcell(5)).peak(indf);

                        %Mean Rate
                        tmp.fmean = f.output{1}(famcell(5)).meanrate(indf);

                        %Activation probability and SWR rate
                        if ~isempty(decodefilter.output{1}(famcell(5)).activationprobability)
                            tmp.activationf = decodefilter.output{1}(famcell(5)).activationprobability(indf);
                            tmp.swrratef = decodefilter.output{1}(famcell(5)).SWRrate(indf);
                            tmp.swrlocf = decodefilter.output{1}(famcell(5)).SWRlocation(indf,[2 3]);
                        else
                            tmp.activationf = NaN;
                            tmp.swrratef = NaN;
                            tmp.swrlocf = nan(1,2);
                        end        
                    else
                        famcell = [];
                        tmp.peakf = nan(1,1); tmp.fmean = nan(1,1); tmp.activationf = nan(1,1); tmp.swrratef = nan(1,1); tmp.swrlocf = nan(1,2);
                    end
                else
                    tmp.peakf = nan(1,1); tmp.fmean = nan(1,1);  tmp.activationf = nan(1,1); tmp.swrratef = nan(1,1); tmp.swrlocf = nan(1,2);
                end
                if ~isempty(famcell) || ~isempty(novcell)
                    out.novel(cellindex) = tmp;
                end
            end
            clear tmp ind indn indf
        end
    end
    %Get rid of novel cells that have already been accounted for
    if ~isempty(nov) && ~isempty(sup)
        for i = 1:size(nov,1)
            if rowfind(nov(i,[1 3 4]),sup(:,[1 3 4]))
                %There is a matching entry in the super novel case, take 
                %out this index from the novel cell only category
                nov(i,:) = NaN;
            end
        end
        nov(isnan(nov(:,1)),:) = [];
    end

    % Go through each remaining novel cell
    for cellindex = 1:size(nov,1)
        novcell = nov(cellindex,:);
        ind = rowfind(novcell(1:4),f.output{2}(novcell(5)).index);
        
        if ~isnan(f.output{2}(novcell(5)).quadrants(ind))
            tmp.index = novcell([1 3 4]);
            
            %Preallocate for supernovel
            tmp.types = NaN;
            tmp.peaks = NaN;
            tmp.smean = NaN;
            tmp.activations = NaN;
            tmp.swrrates = NaN;
            tmp.swrlocs = nan(1,2);
            
            %Determine the quadrant type
            tmp.typen = f.output{2}(novcell(5)).type(ind);
            
            %Look at peak firing rate for novel only cell
            tmp.peakn = f.output{2}(novcell(5)).peak(ind);

            %Determine the mean firing rate
            tmp.nmean = f.output{2}(novcell(5)).meanrate(ind);

            %Activation probability and SWR rate
            if ~isempty(decodefilter.output{2}(novcell(5)).activationprobability)
                tmp.activationn = decodefilter.output{2}(novcell(5)).activationprobability(ind);
                tmp.swrraten = decodefilter.output{2}(novcell(5)).SWRrate(ind);
                tmp.swrlocn = decodefilter.output{2}(novcell(5)).SWRlocation(ind,[2 3]);
            else
                tmp.activationn = NaN;
                tmp.swrraten = NaN;
                tmp.swrlocn = nan(1,2);
            end

            %Find the matching familiar cell
            famcell = fam(find(rowfind(fam(:,1:4),novcell([1 6 3 4]))),:);
            if ~isempty(famcell)
                indf = rowfind(famcell(1:4),f.output{1}(famcell(5)).index);
                if ~isnan(f.output{1}(famcell(5)).peak(indf));
                    
                    %Peak rate
                    tmp.peakf = f.output{1}(famcell(5)).peak(indf);

                    %Mean rate
                    tmp.fmean = f.output{1}(famcell(5)).meanrate(indf);

                    %Activation probability and SWR rate
                    if ~isempty(decodefilter.output{1}(famcell(5)).activationprobability)
                        tmp.activationf = decodefilter.output{1}(famcell(5)).activationprobability(indf);
                        tmp.swrratef = decodefilter.output{1}(famcell(5)).SWRrate(indf);
                        tmp.swrlocf = decodefilter.output{1}(famcell(5)).SWRlocation(indf,[2 3]);
                    else
                        tmp.activationf = NaN;
                        tmp.swrratef = NaN;
                        tmp.swrlocf = nan(1,2);
                    end

                    out.novel(cellindex) = tmp;
                end
            end
            clear tmp ind indf
        end
    end

    
else
    out.novel = [];
end