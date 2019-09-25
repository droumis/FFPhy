function dioprocess(daydirect, animdirect, fileprefix, daynum, varargin)

% function dioprocess(daydirect, animdirect, fileprefix, daynum, varargin)
%
%DAYDIRECT -- folder name where the day's raw data is stored
%ANIMDIRECT -- the path to where the animal's processed data will be stored -- example '/data99/student/Mil'
%FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
%DAYNUM -- the day number for the experiment (starting with 1)

TIMESTAMPRATE = 10000;

[otherArgs] = procOptions(varargin);

times = gettimes(fullfile(daydirect,'times.mat'));

for epoch = 1:length(times)
    if ~isempty(times(epoch).starttime)
        task{daynum}{epoch}.epoch = epoch;
        task{daynum}{epoch}.day = daynum;
    end
end

save(fullfile(animdirect,sprintf('%stask%02d.mat',fileprefix,daynum)),'task');
    
