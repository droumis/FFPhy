function output = calccoactivationstats(g,decodefilter)

nov = []; sup = [];
for an = 1:length(g)
    tmpout = calccoactivationstats2(g(an),decodefilter(an));
    if isfield(tmpout,'novel')
        nov = [nov tmpout.novel];
    end
    if isfield(tmpout,'supernovel')
        sup = [sup tmpout.supernovel];
    end
end


output.corr = []; output.overlap = []; output.peakn = []; 
output.peakf = []; output.type = []; output.nmean = []; output.fmean = []; 
output.locn = []; output.locf = []; output.quad = []; output.preference = [];
output.activation = []; output.swrrate = [];
output.activationf = []; output.swrratef = [];

for i = 1:length(nov)
    if ~isempty(nov(i).nrate)
        %Calculate the correlation between novel and familiar
        tmp = [reshape(nov(i).nrate,size(nov(i).nrate,1)*size(nov(i).nrate,2),1) ...
            reshape(nov(i).frate,size(nov(i).frate,1)*size(nov(i).frate,2),1)];
        invalid = isnan(tmp(:,1)) | isnan(tmp(:,2));
        output.corr = [output.corr; corr(tmp(~invalid,1),tmp(~invalid,2))];
        overlap1 = 2*sum(sum(nov(i).nrate > 1 & nov(i).frate>1));
        overlap2 = (sum(sum(nov(i).nrate>1&nov(i).frate>=0))+sum(sum(nov(i).nrate>=0 & nov(i).frate>1)));
        output.overlap = [output.overlap; overlap1./overlap2];
        output.type = [output.type; nov(i).type];
        output.quad = [output.quad; nov(i).quad];
        output.preference = [output.preference; nov(i).preference];
        output.peakn = [output.peakn; nov(i).peakn];  output.nmean = [output.nmean; nov(i).nmean];
        output.peakf = [output.peakf; nov(i).peakf];  output.fmean = [output.fmean; nov(i).fmean];
        output.locn = [output.locn; nov(i).nlocation];    output.locf = [output.locf; nov(i).flocation];
        output.activation = [output.activation; nov(i).activationn];
        output.activationf = [output.activationf; nov(i).activationf];
        output.swrrate = [output.swrrate; nov(i).swrraten];
        output.swrratef = [output.swrratef; nov(i).swrratef];
    end
end
for i = 1:length(sup)
    if ~isempty(sup(i).srate) && ~isempty(sup(i).frate)
        %Calculate the correlation between super novel and familiar
        tmp = [reshape(sup(i).srate,size(sup(i).srate,1)*size(sup(i).srate,2),1) ...
            reshape(sup(i).frate,size(sup(i).frate,1)*size(sup(i).frate,2),1)];
        invalid = isnan(tmp(:,1)) | isnan(tmp(:,2));
        output.corr = [output.corr; corr(tmp(~invalid,1),tmp(~invalid,2))];
        overlap1 = 2*sum(sum(sup(i).srate > 1 & sup(i).frate>1));
        overlap2 = (sum(sum(sup(i).srate>1&sup(i).frate>=0))+sum(sum(sup(i).srate>=0 & sup(i).frate>1)));
        output.overlap = [output.overlap; overlap1./overlap2];
        output.type = [output.type; sup(i).type];
        output.quad = [output.quad; sup(i).quad];
        output.preference = [output.preference; sup(i).preference];
        output.peakn = [output.peakn; sup(i).peaks];	output.nmean = [output.nmean; sup(i).nmean];
        output.peakf = [output.peakf; sup(i).peakf];	output.fmean = [output.fmean; sup(i).fmean];
        output.locn = [output.locn; sup(i).slocation];    output.locf = [output.locf; sup(i).flocation];
        output.activation = [output.activation; sup(i).activations];
        output.activationf = [output.activationf; sup(i).activationf];
        output.swrrate = [output.swrrate; sup(i).swrrates];
        output.swrratef = [output.swrratef; sup(i).swrratef];
    end
    if ~isempty(sup(i).srate) && ~isempty(sup(i).nrate) && supnovel
        %Calculate the correlation between super novel and novel
        tmp = [reshape(sup(i).srate,size(sup(i).srate,1)*size(sup(i).srate,2),1) ...
            reshape(sup(i).nrate,size(sup(i).nrate,1)*size(sup(i).nrate,2),1)];
        invalid = isnan(tmp(:,1)) | isnan(tmp(:,2));
        output.corr = [output.corr; corr(tmp(~invalid,1),tmp(~invalid,2))];
        overlap1 = 2*sum(sum(sup(i).srate > 1 & sup(i).nrate>1));
        overlap2 = (sum(sum(sup(i).srate>1&sup(i).nrate>=0))+sum(sum(sup(i).srate>=0 & sup(i).nrate>1)));
        output.overlap = [output.overlap; overlap1./overlap2];
        output.type = [output.type; sup(i).type];
        output.quad = [output.quad; sup(i).quad];
        output.preference = [output.preference; sup(i).preference];
        output.peakn = [output.peakn; sup(i).peaks];  output.nmean = [output.nmean; sup(i).nmean];
        output.peakf = [output.peakf; sup(i).peakf];  output.fmean = [output.fmean; sup(i).fmean];
        output.locn = [output.locn; sup(i).slocation];    output.locf = [output.locf; sup(i).nlocation];
        output.activation = [output.activation; sup(i).activations];
        output.activationf = [output.activationf; sup(i).activationn];
        output.swrrate = [output.swrrate; sup(i).swrrates];
        output.swrratef = [output.swrratef; sup(i).swrraten];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = calccoactivationstats2(f,decodefilter)

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

    % For each novel cell, calculate the peak and average difference in firing rate
    for cellindex = 1:size(nov,1)

        %Look at peak firing rate for novel cell
        novcell = nov(cellindex,:);
        ind = rowfind(novcell(1:4),f.output{2}(novcell(5)).index);
        tmp.peakn = f.output{2}(novcell(5)).peak(ind);

        %Determine the quadrant and quadrant type
        tmp.type = f.output{2}(novcell(5)).type(ind);
        tmp.quad = f.output{2}(novcell(5)).quadrants(ind);
        tmp.preference = f.output{2}(novcell(5)).preference;
        
        %Determine the firing rate map
        nx = f.output{2}(novcell(5)).xticks; ny = f.output{2}(novcell(5)).yticks;
        tmp.nrate = f.output{2}(novcell(5)).rate{ind};
        
        %Determine the mean firing rate
        tmp.nmean = f.output{2}(novcell(5)).meanrate(ind);
        
        %Location
        tmp.nlocation = f.output{2}(novcell(5)).location{ind}(f.output{2}(novcell(5)).quadrants(ind),:);
        
        %Activation probability and SWR rate
        if ~isempty(decodefilter.output{2}(novcell(5)).activationprobability)
            tmp.activationn = decodefilter.output{2}(novcell(5)).activationprobability(ind);
            tmp.swrraten = decodefilter.output{2}(novcell(5)).SWRrate(ind);
        else
            tmp.activationn = NaN;
            tmp.swrraten = NaN;
        end
        
        %Find the matching familiar cell
        famcell = fam(find(rowfind(fam(:,1:4),novcell([1 6 3 4]))),:);
        if ~isempty(famcell)
            indf = rowfind(famcell(1:4),f.output{1}(famcell(5)).index);
            if ~all(isnan(f.output{1}(famcell(5)).peak(indf)));
                %Peak rate
                tmp.peakf = f.output{1}(famcell(5)).peak(indf);

                %Determine the firing rate map
                fx = f.output{1}(famcell(5)).xticks; fy = f.output{1}(famcell(5)).yticks;
                tmp.frate = f.output{1}(famcell(5)).rate{indf};

                %Mean rate
                tmp.fmean = f.output{1}(famcell(5)).meanrate(indf);

                %Location
                tmp.flocation = f.output{1}(famcell(5)).location{indf}(f.output{2}(novcell(5)).quadrants(ind),:);

                %Activation probability and SWR rate
                if ~isempty(decodefilter.output{1}(famcell(5)).activationprobability)
                    tmp.activationf = decodefilter.output{1}(famcell(5)).activationprobability(indf);
                    tmp.swrratef = decodefilter.output{1}(famcell(5)).SWRrate(indf);
                else
                    tmp.activationf = NaN;
                    tmp.swrratef = NaN;
                end
                
                if ~isequal(fx,nx)
                    [x a b] = intersect(nx,fx);
                    tmp.nrate = tmp.nrate(a,:); tmp.frate = tmp.frate(b,:);
                end
                if ~isequal(fy,ny)
                    [y a b] = intersect(ny,fy);
                    tmp.nrate = tmp.nrate(:,a); tmp.frate = tmp.frate(:,b);
                end
                out.novel(cellindex) = tmp;
            end
        end
        clear tmp ind indf
    end


    % For each super novel cell, calculate the peak and average difference in firing rate
    for cellindex = 1:size(sup,1)
        
        %Peak Rate
        supcell = sup(cellindex,:);
        ind = rowfind(supcell(1:4),f.output{3}(supcell(5)).index);
        tmp.peaks = f.output{3}(supcell(5)).peak(ind);

        %Determine the quadrant and quadrant type
        tmp.type = f.output{3}(supcell(5)).type(ind);
        tmp.quad = f.output{3}(supcell(5)).quadrants(ind);
        tmp.preference = f.output{3}(supcell(5)).preference;
        
        %Determine the firing rate map
        sx = f.output{3}(supcell(5)).xticks; sy = f.output{3}(supcell(5)).yticks;
        tmp.srate = f.output{3}(supcell(5)).rate{ind};

        %Mean Rate
        tmp.smean = f.output{3}(supcell(5)).meanrate(ind);

        %Location
        tmp.slocation = f.output{3}(supcell(5)).location{ind}(f.output{3}(supcell(5)).quadrants(ind),:);
        
        %Activation probability and SWR rate
        if ~isempty(decodefilter.output{3}(supcell(5)).activationprobability)
            tmp.activations = decodefilter.output{3}(supcell(5)).activationprobability(ind);
            tmp.swrrates = decodefilter.output{3}(supcell(5)).SWRrate(ind);
        else
            tmp.activations = NaN;
            tmp.swrrates = NaN;
        end
               
        %Ignore if the rate is not defined for super novel cell
        if ~isempty(tmp.srate)
            %Find the matching novel cell
            novcell = nov(find(rowfind(nov(:,1:4),supcell([1 7 3 4]))),:);
             
            if ~isempty(novcell)
                indn =rowfind(novcell(1:4),f.output{2}(novcell(5)).index);
                if ~all(isnan(f.output{2}(novcell(5)).peak(indn)));

                    %Peak rate
                    tmp.peakn = f.output{2}(novcell(5)).peak(indn);

                    %Determine the firing rate map
                    nx = f.output{2}(novcell(5)).xticks; ny = f.output{2}(novcell(5)).yticks;
                    tmp.nrate = f.output{2}(novcell(5)).rate{indn};
                    if ~isequal(sx,nx)
                        [x a b] = intersect(nx,sx);
                        tmp.nrate = tmp.nrate(a,:); tmp.srate = tmp.srate(b,:);	sx = x;
                    end
                    if ~isequal(sy,ny)
                        [y a b] = intersect(ny,sy);
                        tmp.nrate = tmp.nrate(:,a); tmp.srate = tmp.srate(:,b);	sy = y;
                    end
                    else
                        novcell = [];
                        tmp.peakn = nan(1,1); tmp.qn = nan(1,1);  tmp.nrate = []; tmp.nmean = nan(1,1);
                    end
                    
                    %Mean rate
                    tmp.nmean = f.output{2}(novcell(5)).meanrate(indn);

                    %Location
                    tmp.nlocation = f.output{2}(novcell(5)).location{indn}(f.output{3}(supcell(5)).quadrants(ind),:);
                    
                    %Activation probability and SWR rate
                    if ~isempty(decodefilter.output{2}(novcell(5)).activationprobability)
                        tmp.activationn = decodefilter.output{2}(novcell(5)).activationprobability(indn);
                        tmp.swrraten = decodefilter.output{2}(novcell(5)).SWRrate(indn);
                    else
                        tmp.activationn = NaN;
                        tmp.swrraten = NaN;
                    end

            else
                tmp.peakn = nan(1,1);  tmp.nrate = []; tmp.nmean = nan(1,1); tmp.nlocation = nan(1,2); tmp.activationn = nan(1,1); tmp.swrraten = nan(1,1);
            end

            %Find the matching familiar cell 
            famcell = fam(find(rowfind(fam(:,1:4),supcell([1 6 3 4]))),:);
            if ~isempty(famcell)
                indf = rowfind(famcell(1:4),f.output{1}(famcell(5)).index);
                if ~all(isnan(f.output{1}(famcell(5)).peak(indf)));
                    %Peak Rate
                    tmp.peakf = f.output{1}(famcell(5)).peak(indf);

                    %Determine the firing rate map
                    fx = f.output{1}(famcell(5)).xticks; fy = f.output{1}(famcell(5)).yticks;
                    tmp.frate = f.output{1}(famcell(5)).rate{indf};

                    if ~isequal(fx,sx)
                        [x a b] = intersect(sx,fx);
                        tmp.srate = tmp.srate(a,:); tmp.frate = tmp.frate(b,:);
                        if ~isempty(tmp.nrate)
                            tmp.nrate = tmp.nrate(a,:);
                        end
                    end
                    if ~isequal(fy,sy)
                        [y a b] = intersect(sy,fy);
                        tmp.srate = tmp.srate(:,a); tmp.frate = tmp.frate(:,b);
                        if ~isempty(tmp.nrate)
                            tmp.nrate = tmp.nrate(:,a);
                        end
                    end
                    
                    %Mean Rate
                    tmp.fmean = f.output{1}(famcell(5)).meanrate(indf);
                    
                    %Location
                    tmp.flocation = f.output{1}(famcell(5)).location{indf}(f.output{3}(supcell(5)).quadrants(ind),:);
                    
                    %Activation probability and SWR rate
                    if ~isempty(decodefilter.output{1}(famcell(5)).activationprobability)
                        tmp.activationf = decodefilter.output{1}(famcell(5)).activationprobability(indf);
                        tmp.swrratef = decodefilter.output{1}(famcell(5)).SWRrate(indf);
                    else
                        tmp.activationf = NaN;
                        tmp.swrratef = NaN;
                    end        
                else
                    famcell = [];
                    tmp.peakf = nan(1,1); tmp.frate = []; tmp.fmean = nan(1,1); tmp.flocation = nan(1,2); tmp.activationf = nan(1,1); tmp.swrratef = nan(1,1);
                end
            else
                tmp.peakf = nan(1,1); tmp.frate = []; tmp.fmean = nan(1,1);  tmp.flocation = nan(1,2); tmp.activationf = nan(1,1); tmp.swrratef = nan(1,1);
            end
            if ~isempty(famcell) || ~isempty(novcell)
                out.supernovel(cellindex) = tmp;
            end
        end
        clear tmp ind indn indf
    end
else
    out.novel = []; out.supernovel = [];
end