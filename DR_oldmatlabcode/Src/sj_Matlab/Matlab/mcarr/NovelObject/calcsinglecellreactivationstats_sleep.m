function output = calcsinglecellreactivationstats_sleep(g,decodefilter)

rates = []; out = [];
for an = 1:length(g)
    tmpout = calcfiringratestats2(g(an),decodefilter(an));
    out = [out tmpout];
end

output.typef = []; output.typen = []; output.types = [];
output.peak = []; output.activationb = []; 
output.activationf = []; output.activationn = []; output.activations = [];

for i = 1:length(out)
    if ~isempty(out{i})
        output.typef = [output.typef; out{i}.type(1)];
        output.typen = [output.typen; out{i}.type(2)];
        if isfield(out{i},'sup')
            output.types = [output.types; out{i}.type(3)];
        else
            output.types = [output.types; NaN];
        end
        output.peak = [output.peak; out{i}.peak];
        
        if isfield(out{i},'baseline')
            output.activationb = [output.activationb; out{i}.baseline];
        else
            output.activationb = [output.activationb; NaN];
        end
        if isfield(out{i},'fam')
            output.activationf = [output.activationf; out{i}.fam];
        else
            output.activationf = [output.activationf; NaN];
        end
          if isfield(out{i},'nov')
            output.activationn = [output.activationn; out{i}.nov];
        else
            output.activationn = [output.activationn; NaN];
        end
          if isfield(out{i},'sup')
            output.activations = [output.activations; out{i}.sup];
        else
            output.activations = [output.activations; NaN];
        end
  

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = calcfiringratestats2(f,d)

%Initialize output
out = [];

% Determine the index for familiar, novel, and super novel sessions
fam = [];   nov = [];   sup = [];
for i = 1:length(f.output{1});
    tmp = f.output{1}(i).index;
    if ~isempty(tmp)
        tmp(:,5) = i;
    end
    fam = [fam; tmp];
end

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

%Create a grand cell index and list of quadrant types
cellindex = [fam zeros(size(fam,1),13)];

%For each cell, come up with its firing rate map
s = {}; o = {}; binx = {}; biny = {};
for c = 1:size(fam,1)
    famcell = fam(c,:);
    ind = rowfind(famcell(1:4),f.output{1}(famcell(5)).index);
    cellindex(c,6) = ind;
     
    if ~isempty(f.output{1}(famcell(5)).spikes{ind})
        s{c} = f.output{1}(famcell(5)).spikes{ind};
    else
        s{c} = zeros(size(f.output{1}(famcell(5)).occupancy));
    end
    o{c} = f.output{1}(famcell(5)).occupancy;
    binx{c} = f.output{1}(famcell(5)).xticks;
    biny{c} = f.output{1}(famcell(5)).yticks;

end

%Add novel cells that have already been accounted for
if ~isempty(nov) && ~isempty(fam)
    for i = 1:size(fam,1)
        if rowfind(fam(i,[1 3 4]),nov(:,[1 3 4]))
            %There is a matching entry in the novel case, add this
            %spike and occupancy information to the cell
            novcell = nov(rowfind(fam(i,[1 3 4]),nov(:,[1 3 4])),:);
            ind = rowfind(novcell(1:4),f.output{2}(novcell(5)).index);
            if ~isempty(f.output{2}(novcell(5)).spikes{ind})
            
                cellindex(i,7:12) = [novcell(1:5) ind];
                 
                if isequal(f.output{2}(novcell(5)).xticks,binx{i}) && isequal(f.output{2}(novcell(5)).yticks,biny{i})
                    s{i} = s{i} + f.output{2}(novcell(5)).spikes{ind};
                    o{i} = o{i} + f.output{2}(novcell(5)).occupancy;
                else
                    x = f.output{2}(novcell(5)).xticks;
                    y = f.output{2}(novcell(5)).yticks;
                    [cx ax bx] = intersect(x,binx{i});
                    [cy ay by] = intersect(y,biny{i});
                    binx{i} = cx;
                    biny{i} = cy;

                    s{i} = s{i}(bx,by) + f.output{2}(novcell(5)).spikes{ind}(ax,ay);
                    o{i} = o{i}(bx,by) + f.output{2}(novcell(5)).occupancy(ax,ay);
                end
            end
        end
    end
end

%Add supernovel cells that have already been accounted for
if ~isempty(sup) && ~isempty(fam)
    for i = 1:size(fam,1)
        if rowfind(fam(i,[1 3 4]),sup(:,[1 3 4]))
            %There is a matching entry in the supernovel case, add this
            %spike and occupancy information to the cell
            supcell = sup(rowfind(fam(i,[1 3 4]),sup(:,[1 3 4])),:);
            ind = rowfind(supcell(1:4),f.output{3}(supcell(5)).index);
            if ~isempty(f.output{3}(supcell(5)).spikes{ind})
                cellindex(i,13:18) = [supcell(1:5) ind];
                
                if isequal(f.output{3}(supcell(5)).xticks,binx{i}) && isequal(f.output{3}(supcell(5)).yticks,biny{i})
                    s{i} = s{i} + f.output{3}(supcell(5)).spikes{ind};
                    o{i} = o{i} + f.output{3}(supcell(5)).occupancy;
                else
                    x = f.output{3}(supcell(5)).xticks;
                    y = f.output{3}(supcell(5)).yticks;
                    [cx ax bx] = intersect(x,binx{i});
                    [cy ay by] = intersect(y,biny{i});
                    binx{i} = cx;
                    biny{i} = cy;

                    s{i} = s{i}(bx,by) + f.output{3}(supcell(5)).spikes{ind}(ax,ay);
                    o{i} = o{i}(bx,by) + f.output{3}(supcell(5)).occupancy(ax,ay);
                end
            end
        end
    end
end

%Add novel cells that have not already been accounted for
if ~isempty(nov)
    for i = 1:size(nov,1)
        if ~rowfind(nov(i,[1 3 4]),cellindex(:,[1 3 4]))
            %There is a cell in the novel case not in familiar, add this
            %spike and occupancy information to the cell
            novcell = nov(i,:);
            ind = rowfind(novcell(1:4),f.output{2}(novcell(5)).index);
            if ~isempty(f.output{2}(novcell(5)).spikes{ind})
                cellindex(end+1,:) = 0;
                cellindex(end,7:12) = [novcell(1:5) ind];
                
                s{end+1} = f.output{2}(novcell(5)).spikes{ind};
                o{end+1} = f.output{2}(novcell(5)).occupancy;
                binx{end+1} = f.output{2}(novcell(5)).xticks;
                biny{end+1} = f.output{2}(novcell(5)).yticks;

                %Add information from supcell if it was recorded
                if ~isempty(sup)
                    ind = rowfind(novcell([1 3 4]),sup(:,[1 3 4]));
                    if ind~=0
                        supcell = sup(rowfind(novcell([1 3 4]),sup(:,[1 3 4])),:);
                        ind = rowfind(supcell(1:4),f.output{3}(supcell(5)).index);
                        cellindex(end,13:18) = [supcell(1:5) ind];
                
                        if isequal(f.output{3}(supcell(5)).xticks,binx{i}) && isequal(f.output{3}(supcell(5)).yticks,biny{i})
                            s{end} = s{end} + f.output{3}(supcell(5)).spikes{ind};
                            o{end} = o{end} + f.output{3}(supcell(5)).occupancy;
                        else
                            x = f.output{3}(supcell(5)).xticks;
                            y = f.output{3}(supcell(5)).yticks;
                            [cx ax bx] = intersect(x,binx{end});
                            [cy ay by] = intersect(y,biny{end});
                            binx{end} = cx;
                            biny{end} = cy;

                            s{end} = s{end}(bx,by) + f.output{3}(supcell(5)).spikes{ind}(ax,ay);
                            o{end} = o{end}(bx,by) + f.output{3}(supcell(5)).occupancy(ax,ay);
                        end
                    end
                end
            end
        end
    end
end

%Add supernovel cells that have not already been accounted for
if ~isempty(sup)
    for i = 1:size(sup,1)
        if ~rowfind(sup(i,[1 3 4]),cellindex(:,[1 3 4]))
            %There is a cell in the novel case not in familiar, add this
            %spike and occupancy information to the cell
            supcell = sup(i,:);
            ind = rowfind(supcell(1:4),f.output{3}(supcell(5)).index);
            if ~isempty(f.output{3}(supcell(5)).spikes{ind})
                cellindex(end+1,:) = 0;
                cellindex(end,13:18) = [supcell(1:5) ind];
                
                s{end+1} = f.output{3}(supcell(5)).spikes{ind};
                o{end+1} = f.output{3}(supcell(5)).occupancy;
                binx{end+1} = f.output{3}(supcell(5)).xticks;
                biny{end+1} = f.output{3}(supcell(5)).yticks;
            end
        end
    end
end

%Now that you have the spikes and occupancy for each cell, construct a rate
%map for each cell and describe what quadrant it is in
g = gaussian2(8,24);
rate = cell(size(s)); cellinfo = cell(size(s));
linearcoords = []; i = 1;
while isempty(linearcoords)
    linearcoords = f.output{1}(i).linearcoords;
    if isempty(linearcoords)
        i = i+1;
    end
end

for c = 1:size(cellindex,1)
    
    %Make sure you have enough spikes to say something
    if sum(sum(s{c})) > 50
        rate{c} = conv2(s{c},g,'same')./conv2(o{c},g,'same');
        
        %Determine whether the cell is a place cell
        if max(max(rate{c})) >= 1
            
            %Determine which quadrant the place cell is in
            [maxy yind] = max(rate{c});
            [maxxy xind] = max(maxy);
            peak = [binx{c}(yind(xind)) biny{c}(xind)];
            
            if cellindex(c,5) ~=0
                quadmapping = f.output{1}(cellindex(c,5)).quadrants;
                preference = f.output{1}(cellindex(c,5)).preference;
            else
                quadmapping = nan(4,1);
                preference = NaN;
            end
            if cellindex(c,11) ~=0
                quadmapping = [quadmapping f.output{2}(cellindex(c,11)).quadrants];
                preference = [preference f.output{2}(cellindex(c,11)).preference];
            else
                quadmapping = [quadmapping nan(4,1)];
                preference = [preference NaN];
            end
            if cellindex(c,17) ~=0
                quadmapping = [quadmapping f.output{3}(cellindex(c,17)).quadrants];
                preference = [preference f.output{3}(cellindex(c,17)).preference];
            else
                quadmapping = [quadmapping nan(4,1)];
                preference = [preference NaN];
            end
            
            %Determine whether there was enough data near the peak
            occ = conv2(o{c},g,'same');
            if occ(yind(xind),xind) > 0.001
                if peak(1) <= linearcoords(1) && peak(2) >= linearcoords(2)
                    q = 1;
                elseif peak(1) > linearcoords(1) && peak(2) >= linearcoords(2)
                    q = 2;
                elseif peak(1) > linearcoords(1) && peak(2) < linearcoords(2)
                    q = 3;
                elseif peak(1) <= linearcoords(1) && peak(2) < linearcoords(2)
                    q = 4;
                end
                cellinfo{c}.quadrant = q;
                cellinfo{c}.peak = maxxy;
                cellinfo{c}.type = quadmapping(q,:);
                cellinfo{c}.index = cellindex(c,:);
                cellinfo{c}.rate = rate{c};
                cellinfo{c}.preference = preference;
            else
                rate{c} = [];
            end
        else
            rate{c} = [];
        end
    end
end

clear cellindex rate binx biny o s linearcoords quadmapping g peak maxy yind maxxy xind occ

%Compute single cell reactivation

% Determine the index for baseline, familiar, novel, and super novel sessions
bas = []; fam = [];   nov = [];   sup = [];
for i = 1:length(d.output{1});
    tmp = d.output{1}(i).index; tmp(:,5) = i; tmp(:,6) = 1:size(tmp,1);
    bas = [bas; tmp];
end

for i = 1:length(d.output{2});
    tmp = d.output{2}(i).index; tmp(:,5) = i; tmp(:,6) = 1:size(tmp,1);
    fam = [fam; tmp];
end

for i = 1:length(d.output{3});
    tmp = d.output{3}(i).index; tmp(:,5) = i; tmp(:,6) = 1:size(tmp,1);
    nov = [nov; tmp];
end

if length(d.output) == 4
    for i = 1:length(d.output{4});
        tmp = d.output{4}(i).index; tmp(:,5) = i; tmp(:,6) = 1:size(tmp,1);
        sup = [sup; tmp];
    end
end
clear tmp

%Go through each place cell and determine activation probability
for c = 1:length(cellinfo)
    if ~isempty(cellinfo{c})
        %Determine which baseline cell corresponds to place cell
        ind =rowfind(cellinfo{c}.index([1 3:4]),bas(:,[1 3:4]));
        if ind~=0 && ~isempty(d.output{1}(bas(ind,5)).index)
            cellinfo{c}.baseline = d.output{1}(bas(ind,5)).activationprobability(bas(ind,6));
        end    
        %Determine which fam cell corresponds to place cell
        ind =rowfind(cellinfo{c}.index([1 3:4]),fam(:,[1 3:4]));
        if ind ~=0 && ~isempty(d.output{2}(fam(ind,5)).index)
            cellinfo{c}.fam = d.output{2}(fam(ind,5)).activationprobability(fam(ind,6));
        end    
        %Determine which nov cell corresponds to place cell
        ind =rowfind(cellinfo{c}.index([1 3:4]),nov(:,[1 3:4]));
        if ind~=0   && ~isempty(d.output{3}(nov(ind,5)).index)
            cellinfo{c}.nov = d.output{3}(nov(ind,5)).activationprobability(nov(ind,6));
        end
        %Determine which sup cell corresponds to place cell
        if ~isempty(sup)
            ind =rowfind(cellinfo{c}.index([1 3:4]),sup(:,[1 3:4]));
            if ind~=0 && ~isempty(d.output{4}(sup(ind,5)).index)
                cellinfo{c}.sup = d.output{4}(sup(ind,5)).activationprobability(sup(ind,6));
            end
        end
    end
end

out = cellinfo;