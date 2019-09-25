function createripstruct(directoryname, fileprefix, days, varargin)
% function createripstruct(directoryname, fileprefix, days, varargin)
%   This function goes through each day - epoch and creates the rip, rips,
%   and ripc structures. These structures contains the information about
%   each ripple and are saved in the animal's processed data directory

%   rip struct:
%   -starttime
%   -endtime
%   -spiketimes
%   -cells that fired during the ripple
%   -time relative to ripple onset for spectrum and coherence measures
%   -average CA1 gamma power
%   -average CA1 ripple power
%   -average CA3 gamma power
%   -average CA1-CA1 gamma coherence
%   -average CA1-CA3 gamma cohernece
%   -average CA3-CA3 gamma coherence

%   rips struct:
%   -starttime
%   -time relative to ripple onset for spectrum and coherence measures
%   -frequency
%   -average CA1 spectrum
%   -average CA3 spectrum

%   ripc struct:
%   -starttime
%   -time relative to ripple onset for spectrum and coherence measures
%   -frequency
%   -average CA1-CA1 coherence
%   -CA1-CA1 phase distribution
%   -average CA1-CA3 cohernece
%   -CA1-CA3 phase distribution
%   -average CA3-CA3 coherence
%   -CA3-CA3 phase distribution   

% set options
savedirectoryname = directoryname;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'bin'
                bin = varargin{option+1};
            case 'savedirectory'
                savedirectoryname = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% check to see if the directory has a trailing '/'
if (directoryname(end) ~= '/')
    warning('directoryname should end with a ''/'', appending one and continuing');
    directoryname(end+1) = '/';
end

days = days(:)';

%Load cellinfo and tetinfo
load(sprintf([directoryname,fileprefix,'cellinfo.mat']))
load(sprintf([directoryname,fileprefix,'tetinfo.mat']))

for day = days
    rip = cell(1,day);
    rips = cell(1,day);
    ripc = cell(1,day);
    
    %Load ripples and spikes for this day
    load(sprintf('%s/%sspikes%02d.mat',directoryname,fileprefix,day))
    load(sprintf('%s/%sripples%02d.mat',directoryname,fileprefix,day))
    
    rip{day} = cell(1,min(length(spikes{day}),7));
    rips{day} = cell(1,min(length(spikes{day}),7));
    ripc{day} = cell(1,min(length(spikes{day}),7));
    
    for epoch = 1:min(length(spikes{day}),7)
        if ~isempty(ripples{day}{epoch})
            %Determine when ripples occur and add starttime and endtime to
            %the rip structure and starttime to the rips and ripc structure
            riptimes = getripples([day epoch], ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1''))','minstd',3);
            rip{day}{epoch}.starttime = riptimes(:,1);
            rip{day}{epoch}.endtime = riptimes(:,2);
            
            rips{day}{epoch}.starttime = riptimes(:,1);
            ripc{day}{epoch}.starttime = riptimes(:,1);
            
            %Determine which cells fire during each ripple and add cell
            %identity and spiketimes to the rip structure
            indices = evaluatefilter(cellinfo{day}{epoch},'isequal($area,''CA1'') | isequal($area,''CA3'') & $meanrate < 7');

            celldata = [];
            for cellcount = 1:size(indices,1)
                index = [day epoch indices(cellcount,:)];
                if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                else
                    spiketimes = [];
                end

                %Find valid spikes
                spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
                if ~isempty(spiketimes)
                    validspikes = find(spikebins);
                    spiketimes = spiketimes(validspikes);
                    spikebins = spikebins(validspikes);
                end

                tmpcelldata = [spiketimes spikebins];
                tmpcelldata(:,3) = cellcount;
                celldata = [celldata; tmpcelldata];
            end
            celldata = sortrows(celldata,1);

            for ripcount = 1:size(riptimes,1)
                rip{day}{epoch}.spiketimes{ripcount} = celldata(celldata(:,2)==ripcount,1);
                rip{day}{epoch}.celldata{ripcount} = indices(celldata(celldata(:,2)==ripcount,3),:);
            end
            clear celldata tmpcelldata spikebins spiketimes indices riptimes

            %Go through all tetrodes with more than one cell and compute
            %average CA1 spectrum, average CA3 spectrum, avereage CA1-CA1
            %coherence, average CA1-CA3 coherence, and average CA3-CA3
            %coherence, and the gamma power and coherence for each ripple
            ca1tetrodes = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''CA1'') & $numcells > 1');
            ca3tetrodes = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''CA3'') & $numcells > 1');

            rip{day}{epoch}.time = [];
            tmp = zeros(91,59,size(rip{day}{epoch}.starttime,1)); count = 0;
            for i = ca1tetrodes'
                load(sprintf('%sEEGnonreference/%sspectrum%02d-%d-%02d.mat',directoryname,fileprefix,day,epoch,i))
                ripind = lookup(rip{day}{epoch}.starttime,spectrum{day}{epoch}{i}.ripples);
                tmp = tmp + spectrum{day}{epoch}{i}.spectrum(:,:,ripind);
                count = count+1;
                if isempty(rip{day}{epoch}.time)
                    rip{day}{epoch}.time = spectrum{day}{epoch}{i}.time;
                    
                    rips{day}{epoch}.time = spectrum{day}{epoch}{i}.time;
                    rips{day}{epoch}.frequency = spectrum{day}{epoch}{i}.frequency;
                    
                    ripc{day}{epoch}.time = spectrum{day}{epoch}{i}.time;
                    ripc{day}{epoch}.frequency = spectrum{day}{epoch}{i}.frequency;
                end
                clear spectrum
            end
            rips{day}{epoch}.ca1_spectrum = tmp./count;
 
            tmp = zeros(91,59,size(rip{day}{epoch}.starttime,1)); count = 0;
            for i = ca3tetrodes'
                load(sprintf('%sEEGnonreference/%sspectrum%02d-%d-%02d.mat',directoryname,fileprefix,day,epoch,i))
                ripind = lookup(rip{day}{epoch}.starttime,spectrum{day}{epoch}{i}.ripples);
                tmp = tmp + spectrum{day}{epoch}{i}.spectrum(:,:,ripind);
                count = count+1;
                if isempty(rip{day}{epoch}.time)
                    rip{day}{epoch}.time = spectrum{day}{epoch}{i}.time;
                    
                    rips{day}{epoch}.time = spectrum{day}{epoch}{i}.time;
                    rips{day}{epoch}.frequency = spectrum{day}{epoch}{i}.frequency;
                    
                    ripc{day}{epoch}.time = spectrum{day}{epoch}{i}.time;
                    ripc{day}{epoch}.frequency = spectrum{day}{epoch}{i}.frequency;
                end
                clear spectrum
            end
            rips{day}{epoch}.ca3_spectrum = tmp./count;
            clear tmp count
            
            %Determine the average CA1 gamma power and ripple power
            gamma_index = lookup([20 50],rips{day}{epoch}.frequency);
            ripple_index = lookup([150 250],rips{day}{epoch}.frequency);
            rip{day}{epoch}.ca1_gamma_power = squeeze(mean(rips{day}{epoch}.ca1_spectrum(:,gamma_index(1):gamma_index(2),:),2));
            rip{day}{epoch}.ca1_ripple_power = squeeze(mean(rips{day}{epoch}.ca1_spectrum(:,ripple_index(1):ripple_index(2),:),2));
            rip{day}{epoch}.ca3_gamma_power = squeeze(mean(rips{day}{epoch}.ca3_spectrum(:,gamma_index(1):gamma_index(2),:),2));

            %Initialize variables
            tmp1 = zeros(length(rip{day}{epoch}.time),59,size(rip{day}{epoch}.starttime,1)); count1 = 0;
            p1 = []; 

            tmp13 = zeros(length(rip{day}{epoch}.time),59,size(rip{day}{epoch}.starttime,1)); count13 = 0;
            p13 = [];

            tmp3 = zeros(length(rip{day}{epoch}.time),59,size(rip{day}{epoch}.starttime,1)); count3 = 0;
            p3 = [];

            %Compute coherence for CA1-CA1 and CA1-CA3 tetrode pairs
            for i = ca1tetrodes'
                load(sprintf('%sEEGnonreference/%scoherence%02d-%d-%02d.mat',directoryname,fileprefix,day,epoch,i))
                for j = ca1tetrodes'
                    if j > i
                        ripind = lookup(rip{day}{epoch}.starttime,coherence{day}{epoch}{i}{j}.ripples);
                        tmp1 = tmp1 + im2double(coherence{day}{epoch}{i}{j}.coherence(:,:,ripind));
                        count1 = count1+1;
                        if isempty(p1)
                            p1 = squeeze(mean(coherence{day}{epoch}{i}{j}.phase(:,2:8,ripind),2));
                        else
                            p1 = cat(3,p1,squeeze(mean(coherence{day}{epoch}{i}{j}.phase(:,2:8,ripind),2)));
                        end
                    end
                end
                for j = ca3tetrodes'
                    ripind = lookup(rip{day}{epoch}.starttime,coherence{day}{epoch}{i}{j}.ripples);
                    tmp13 = tmp13 + im2double(coherence{day}{epoch}{i}{j}.coherence(:,:,ripind));
                    count13 = count13+1;
                    if isempty(p13)
                        p13 = squeeze(mean(coherence{day}{epoch}{i}{j}.phase(:,2:8,ripind),2));
                    else
                    	p13 = cat(3,p13,squeeze(mean(coherence{day}{epoch}{i}{j}.phase(:,2:8,ripind),2)));
                    end
                end
            end
            ripc{day}{epoch}.ca1_ca1_coherence = tmp1./count1;
            ripc{day}{epoch}.ca1_ca1_phase = squeeze(mean(p1,3));
            ripc{day}{epoch}.ca1_ca3_coherence = tmp13./count13;
            ripc{day}{epoch}.ca1_ca3_phase = squeeze(mean(p13,3));
            clear tmp1 p1 count1 coherence tmp13 p13 count13
        
            %Compute coherence for CA3-CA3 tetrode pairs
            for i = ca3tetrodes'
                load(sprintf('%sEEGnonreference/%scoherence%02d-%d-%02d.mat',directoryname,fileprefix,day,epoch,i))
                for j = ca3tetrodes'
                    if j > i
                        ripind = lookup(rip{day}{epoch}.starttime,coherence{day}{epoch}{i}{j}.ripples);
                        tmp3 = tmp3 + im2double(coherence{day}{epoch}{i}{j}.coherence(:,:,ripind));
                        count3 = count3+1;
                        if isempty(p3)
                            p3 = squeeze(mean(coherence{day}{epoch}{i}{j}.phase(:,2:8,ripind),2));
                        else
                        	p3 = cat(3,p3,squeeze(mean(coherence{day}{epoch}{i}{j}.phase(:,2:8,ripind),2)));
                        end
                    end
                end
            end
            ripc{day}{epoch}.ca3_ca3_coherence = tmp3./count3;
            ripc{day}{epoch}.ca3_ca3_phase = squeeze(mean(p3,3));
            clear tmp3 p3 count3 coherence ca1tetrodes ca3tetrodes
            
            %Average CA1-CA1, CA1-CA3, and CA3-CA3 gamma coherence
            rip{day}{epoch}.ca1_ca1_gamma_coherence = squeeze(mean(ripc{day}{epoch}.ca1_ca1_coherence(:,gamma_index(1):gamma_index(2),:),2));
            rip{day}{epoch}.ca1_ca3_gamma_coherence = squeeze(mean(ripc{day}{epoch}.ca1_ca3_coherence(:,gamma_index(1):gamma_index(2),:),2));
            rip{day}{epoch}.ca3_ca3_gamma_coherence = squeeze(mean(ripc{day}{epoch}.ca3_ca3_coherence(:,gamma_index(1):gamma_index(2),:),2));

        else
            rip{day}{epoch} = [];
            rip{day}{epoch}.starttime = [];
            rip{day}{epoch}.endtime = [];
            rip{day}{epoch}.spiketimes = {};
            rip{day}{epoch}.celldata = {};
            rip{day}{epoch}.time = [];
            rip{day}{epoch}.ca1_gamma_power = [];
            rip{day}{epoch}.ca3_gamma_power = [];
            rip{day}{epoch}.ca1_ca1_gamma_coherence = [];
            rip{day}{epoch}.ca1_ca3_gamma_coherence = [];
            rip{day}{epoch}.ca3_ca3_gamma_coherence = [];
            
            rips{day}{epoch} = [];
            rips{day}{epoch}.starttime = [];
            rips{day}{epoch}.time = [];
            rips{day}{epoch}.ca1_spectrum= [];
            rips{day}{epoch}.ca3_spectrum= [];
            
            ripc{day}{epoch} = [];
            ripc{day}{epoch}.starttime = [];
            ripc{day}{epoch}.time = [];
            rip{day}{epoch}.ca1_ca1_coherence= [];
            rip{day}{epoch}.ca1_ca1_phase= [];
            rip{day}{epoch}.ca1_ca3_coherence= [];
            rip{day}{epoch}.ca1_ca3_phase= [];
            rip{day}{epoch}.ca3_ca3_coherence= [];
            rip{day}{epoch}.ca3_ca3_phase= [];
        end
    end
    %save the resulting file
    filename = sprintf('%s/%srip%02d.mat',savedirectoryname,fileprefix,day);
    save(filename,'rip');
    filename = sprintf('%s/%srips%02d.mat',savedirectoryname,fileprefix,day);
    save(filename,'rips');
    filename = sprintf('%s/%sripc%02d.mat',savedirectoryname,fileprefix,day);
    save(filename,'ripc');
end


