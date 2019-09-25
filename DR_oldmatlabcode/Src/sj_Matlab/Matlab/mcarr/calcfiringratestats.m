function output = calcfiringratestats(g)

nov = []; sup = [];
for an = 1:length(g)
    if ~isempty(g(an).output)
        tmpout = calculate(g(an));
        if isfield(tmpout,'novel')
            nov = [nov tmpout.novel];
        end
        if isfield(tmpout,'supernovel')
            sup = [sup tmpout.supernovel];
        end
    end
end


output.corrn = []; output.corrs = [];
output.peakn = []; output.peaks = []; 
output.nmean = []; output.smean = [];
output.typen = []; output.types = [];
output.quadn = []; output.quads = [];
output.npreference = []; output.spreference = [];
output.nlocation = []; output.slocation = [];
for i = 1:length(nov)
    if ~isempty(nov(i).nrate)
        %Calculate the correlation between novel and familiar
        tmp = [reshape(nov(i).nrate,size(nov(i).nrate,1)*size(nov(i).nrate,2),1) ...
            reshape(nov(i).frate,size(nov(i).frate,1)*size(nov(i).frate,2),1)];
        invalid = isnan(tmp(:,1)) | isnan(tmp(:,2));
        output.corrn = [output.corrn; corr(tmp(~invalid,1),tmp(~invalid,2))];
        output.typen = [output.typen; nov(i).type];
        output.quadn = [output.quadn; nov(i).quad];
        output.npreference = [output.npreference; nov(i).preference];
        output.peakn = [output.peakn; nov(i).peakn./nov(i).peakf];  output.nmean = [output.nmean; nov(i).nmean./nov(i).fmean];
        output.nlocation = [output.nlocation(:); nov(i).nlocation(:)];
    end
end
for i = 1:length(sup)
    if ~isempty(sup(i).srate) && ~isempty(sup(i).frate)
        %Calculate the correlation between super novel and familiar
        tmp = [reshape(sup(i).srate,size(sup(i).srate,1)*size(sup(i).srate,2),1) ...
            reshape(sup(i).frate,size(sup(i).frate,1)*size(sup(i).frate,2),1)];
        invalid = isnan(tmp(:,1)) | isnan(tmp(:,2));
        output.corrs = [output.corrs; corr(tmp(~invalid,1),tmp(~invalid,2))];
        output.types = [output.types; sup(i).type];
        output.quads = [output.quads; sup(i).quad];
        output.spreference = [output.spreference; sup(i).preference];
        output.peaks = [output.peaks; sup(i).peaks./sup(i).peakf];	output.smean = [output.smean; sup(i).nmean./sup(i).fmean];
        output.slocation = [output.slocation(:); sup(i).slocation(:)];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = calculate(f)

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
%        tmp.nlocation = f.output{2}(novcell(5)).location{ind}(f.output{2}(novcell(5)).quadrants(ind),:);
        nlocation = f.output{2}(novcell(5)).location{ind}(find(f.output{2}(novcell(5)).quadmapping==1),:);
        
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
                nlocation = median([nlocation; f.output{1}(famcell(5)).location{indf}(find(f.output{2}(novcell(5)).quadmapping==1),:)]);
                tmp.nlocation = sqrt( (nlocation(1)^2) + nlocation(2)^2);
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
        slocation = f.output{3}(supcell(5)).location{ind}(find(f.output{3}(supcell(5)).quadmapping==1),:);
        
        %Ignore if the rate is not defined for super novel cell
        if ~isempty(tmp.srate)
            if ~isempty(nov)
                %Find the matching novel cell
                novcell = nov(find(rowfind(nov(:,1:4),supcell([1 7 3 4]))),:);
            else
                novcell = [];
            end
            
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
%                    tmp.nlocation = f.output{2}(novcell(5)).location{indn}(f.output{3}(supcell(5)).quadrants(ind),:);
                    %tmp.nlocation = f.output{2}(novcell(5)).location{indn}(find(f.output{3}(supcell(5)).quadmapping==1),:);

            else
                tmp.peakn = nan(1,1);  tmp.nrate = []; tmp.nmean = nan(1,1); %tmp.nlocation = nan(1,2);
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
                    slocation = median([slocation; f.output{1}(famcell(5)).location{indf}(find(f.output{3}(supcell(5)).quadmapping==1),:)]);
                    tmp.slocation = sqrt( (slocation(1)^2) + slocation(2)^2);
                else
                    famcell = [];
                    tmp.peakf = nan(1,1); tmp.frate = []; tmp.fmean = nan(1,1);
                end
            else
                tmp.peakf = nan(1,1); tmp.frate = []; tmp.fmean = nan(1,1); tmp.slocation = NaN;
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