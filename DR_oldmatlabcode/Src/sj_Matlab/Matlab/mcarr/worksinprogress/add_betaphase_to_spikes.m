function add_betaphase_to_spikes(directoryname,fileprefix,days, varargin)

%
%Goes through each cell and adds the beta phase recorded on a 
%representative CA3 tetrode to the spikes structure
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
            %load the reference beta file
            ref_tet = evaluatefilter(tetinfo{day}{epoch},ref_tetrode);
            if isempty(ref_tet)
                ref_tet = evaluatefilter(tetinfo{day}{epoch},'isequal($area, ''CA3'')');
                ref_tet = ref_tet(1);
            end

            loadfile = sprintf('%s/EEGnonreference/%seeg%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,ref_tet);
            load(loadfile);
            
            %Filter for beta-gamma
            e = eeg{day}{epoch}{ref_tet}.data;
            rtime = geteegtimes(eeg{day}{epoch}{ref_tet});
            clear eeg
            filterstring = '/home/mcarr/Src/Matlab/Filters/betafilter.mat';
    
            eval(['load ', filterstring]);
            beta = hilbert(filtfilt(betagammafilter,1,e));
            rphase = angle(beta);
            clear e beta filterstring
              
            %Go through each cell and assign the global beta phase to each spike
            cells = evaluatefilter(cellinfo{day}{epoch},'isequal($area,''CA1'') | isequal($area,''CA3'') & $numspikes > 0');
            tets = unique(cells(:,1));

            for t = tets'
                for c = cells(cells(:,1)==t,2)'
                    if ~isempty(spikes{day}{epoch}{t}{c}.data)
                        s = spikes{day}{epoch}{t}{c}.data(:,1);
                        ind = lookup(s,rtime);
                        spikes{day}{epoch}{t}{c}.beta_phase = rphase(ind);
                    end
                end
            end
        end
    end
    
    %Save spikes file
    savefile = sprintf('%s/%sspikes%02d.mat', directoryname, fileprefix,day);
    save(savefile,'spikes');
        
end