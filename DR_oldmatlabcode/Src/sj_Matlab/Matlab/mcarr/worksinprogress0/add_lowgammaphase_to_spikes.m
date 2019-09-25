function add_lowgammaphase_to_spikes(directoryname,fileprefix,days, varargin)

%
%Goes through each cell and adds both the local gamma phase and the current
%gamma phase recorded on a representative CA3 tetrode to the spikes
%structure
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 
%
%options -
%
%   ref_tetrode: a filter command describing which tetrode to use for the
%       global gamma phase. Default: 'isequal($area,''CA3'') & $representative = 1'

ref_tetrode = 'isequal($area,''CA3'') & $maxcell == 1';

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'ref_tetrode'
            ref_tetrode = varargin{option+1};
    end
end


% check to see if the directory has a trailing '/'
if (directoryname(end) ~= '/')
    warning('directoryname should end with a ''/'', appending one and continuing');
    directoryname(end+1) = '/';
end

days = days(:)';

for day = days
    %load up the spike file
    spikes = loaddatastruct(directoryname,fileprefix, 'spikes',day);
    
    %load the tetinfo and cellinfo file
    tetinfo = loaddatastruct(directoryname,fileprefix,'tetinfo');
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo');
    
    for epoch = 1:min(7,length(spikes{day}))
        if ~isempty(spikes{day}{epoch})
            %load the reference low gamma file
            ref_tet = evaluatefilter(tetinfo{day}{epoch},ref_tetrode);
            if isempty(ref_tet)
                ref_tet = evaluatefilter(tetinfo{day}{epoch},'isequal($area, ''CA3'')');
                ref_tet = ref_tet(1);
            end

            loadfile = sprintf('%s/EEGnonreference/%slowgamma%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,ref_tet);
            load(loadfile);
            rtime = geteegtimes(lowgamma{day}{epoch}{ref_tet});
            rphase = double(lowgamma{day}{epoch}{ref_tet}.data(:,2))./10000;
            clear lowgamma loadfile
            [maxtab mintab] = peakdet(rphase,0.01); clear mintab
            maxtab = maxtab(maxtab(:,2)>pi-0.5,1);
            cum_rphase = zeros(size(rphase));
            if min(diff(rphase(1:maxtab(1,1))))>0
                cum_rphase(1:maxtab(1,1)) = rphase(1:maxtab(1,1));

            else
                cum_rphase(1:maxtab(1,1)) = interp1([1 maxtab(1,1)], [rphase(1) rphase(maxtab(1,1))], 1:maxtab(1,1));
            end
            count = 1;
            for i = 1:length(maxtab)
                if i < length(maxtab)
                    cum_rphase(maxtab(i)+1:maxtab(i+1)) = count*2*pi + interp1([maxtab(i)+1 maxtab(i+1)],[-pi pi],maxtab(i)+1:maxtab(i+1));
                    count = count+1;
                elseif maxtab(i) >= length(rphase)-1
                    cum_rphase(maxtab(i)+1:end) = count*2*pi + rphase(maxtab(i)+1:end);
                else
                    cum_rphase(maxtab(i)+1:end) = count*2*pi + interp1([maxtab(i)+1 length(rphase)],[-pi rphase(end)],maxtab(i)+1:length(rphase));
                end
            end
            
            if min(diff(cum_rphase)) < 0
                slips = find(diff([0; cum_rphase])<0);
                for i = 1:length(slips)
                    if slips(i)-1 > 0 && slips(i)+1 < length(cum_rphase)
                        tmp = interp1([slips(i)-1 slips(i)+1],[cum_rphase(slips(i)-1) cum_rphase(slips(i)+1)],slips(i)-1:slips(i)+1);
                        cum_rphase(slips(i)) = tmp(2);
                    elseif slips(i)-1 == 0
                        cum_rphase(slips(i)) = cum_rphase(slips(i)+1)-0.0001;
                    else
                        cum_rphase(slips(i)) = cum_rphase(slips(i)-1) + 0.0001;
                    end                        
                end
            end
            if min(diff(cum_rphase))<0
                keyboard
            end
            %Go through each tetrode that has cells on it, load the local gamma
            %file, and for each cell, assign both the local and global gamma
            %phase to each spike
            cells = evaluatefilter(cellinfo{day}{epoch},'isequal($area,''CA1'') | isequal($area,''CA3'') & $numspikes > 0');
            tets = unique(cells(:,1));

            for t = tets'
                %Load the low gamma file
                loadfile = sprintf('%s/EEGnonreference/%slowgamma%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,t);
                load(loadfile)
                ltime = geteegtimes(lowgamma{day}{epoch}{t});
                lphase = double(lowgamma{day}{epoch}{t}.data(:,2))./10000;
                clear lowgamma load file

                %Go through each cell on the tetrode and assign local and
                %global gamma phase
                for c = cells(cells(:,1)==t,2)'
                    if ~isempty(spikes{day}{epoch}{t}{c}.data)
                        s = spikes{day}{epoch}{t}{c}.data(:,1);
                        ind = lookup(s,rtime);
                        spikes{day}{epoch}{t}{c}.globalgammaphase = rphase(ind);
                        spikes{day}{epoch}{t}{c}.cumulative_phase = cum_rphase(ind);
                        ind = lookup(s,ltime);
                        spikes{day}{epoch}{t}{c}.lowgammaphase = lphase(ind);
                    end
                end
            end
        end
    end
    
    %Save spikes file
    savefile = sprintf('%s/%sspikes%02d.mat', directoryname, fileprefix,day);
    save(savefile,'spikes');
        
end