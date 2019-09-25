function out = calcnonreferenceeeg(directoryname, fileprefix,days,varargin)
% Adds the reference trace to the EEG file and saves in the EEGnonreference
% folder.
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 

savedirectoryname = directoryname;
%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'savedirectory'
            savedirectoryname = varargin{option+1};
        otherwise
                error(['Option ',varargin{option},' unknown.']);    
    end
end


% check to see if the directory has a trailing '/'
if (directoryname(end) ~= '/')
    warning('directoryname should end with a ''/'', appending one and continuing');
    directoryname(end+1) = '/';
end

days = days(:)';

for day = days
    %Create the list of files that we should add the reference to
    load(sprintf('%s/%stetinfo.mat',directoryname,fileprefix));
    
    tetlist = evaluatefilter(tetinfo{day},'~isequal($area,''Reference'')');
    tetlist = unique(tetlist(:,2));
    tmpflist = [];
    for t = 1:length(tetlist)
        tmp = dir(sprintf('%s/EEG/*eeg%02d-*-%02d.mat',directoryname,day,tetlist(t)));
        tmpflist = [tmpflist; tmp];
    end
    flist = {};
    for i = 1:length(tmpflist)
        flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
    end

    %Go through each file in flist and add the reference to it
    for fnum = 1:length(flist)
        dash = find(flist{fnum} == '-');
        epoch = str2num(flist{fnum}((dash(1)+1):(dash(2)-1)));
        tet = str2num(flist{fnum}((dash(2)+1):(dash(2)+3)));
	
        %Only do this for sleeps (have already done the runs)
        if mod(epoch,2)
            % Load the reference tetrode
            reference_tet = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''Reference'')');
            reflist = dir(sprintf('%s/EEG/*eeg%02d-%d-%02d.mat',directoryname,day,epoch,reference_tet));

            load(sprintf('%s/EEG/%s',directoryname,reflist.name))
            ref = eeg{day}{epoch}{reference_tet};

            load(flist{fnum});
            e = eeg{day}{epoch}{tet}; clear eeg

            eegtimes = geteegtimes(e);
            reftimes = geteegtimes(ref);

            if length(eegtimes)>length(reftimes)
                temp = lookup(reftimes,eegtimes);
                e.data = e.data(temp);
                ref.data = ref.data;
            elseif length(reftimes)>length(eegtimes)
                temp = lookup(eegtimes,reftimes);
                e.data = e.data;
                ref.data = ref.data(temp);
            elseif length(eegtimes)==length(reftimes)
                e.data = e.data;
                ref.data = ref.data;
            end

            e.data = e.data + ref.data;
            eeg{day}{epoch}{tet} = e;
            clear e ref
            
            %save the resulting file
            tmpfilename = sprintf('%sEEGnonreference/%seeg%02d-%d-%02d.mat',...
                savedirectoryname,fileprefix,day,epoch,tet);
            save(tmpfilename,'eeg');
            eval(['clear ',fileprefix,'eeg']);

        end
    end
end