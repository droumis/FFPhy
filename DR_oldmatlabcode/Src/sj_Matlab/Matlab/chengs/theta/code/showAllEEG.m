function showAllEEG
global fmaux
% find all files, extract day, epoch and tetrode
f= dir([fmaux.EEGdir '/*theta*']);
nfiles=size(f,1);
files=cell(nfiles,1);
[files{:}]= deal(f.name);
d=zeros(nfiles,1);
e=zeros(nfiles,1);
t=zeros(nfiles,1);
for i=1:nfiles
    [num c]=sscanf(files{i},[fmaux.prefix 'theta%d-%d-%d']);
    if c~=3 
        error(['EEG filename "' files{i} '" does not conform']);
    end
    d(i)= num(1); e(i)= num(2); t(i)= num(3);
end

oldday= -1;
oldepoch= -1;
% show all 
for i=1:nfiles
%    keyboard
    if(oldepoch > 0 &  e(i)==oldepoch & d(i)==oldday); continue; end
    load([fmaux.EEGdir '/' files{i}]);
    eval(['eeg = ' fmaux.prefix 'theta{d(i)}{e(i)}{t(i)};']);
    nsamp= size(eeg.data,2);
    time=[1:nsamp]/eeg.samprate + eeg.starttime;

    figure(1); clf
    subplot(4,1,1)
    auxplot(time,eeg.data,1,nsamp);
    title(files{i});

    subplot(4,1,2)
    auxplot(time,eeg.data,1,400);
    title('beginning');

    subplot(4,1,3)
    auxplot(time,eeg.data,floor(nsamp/2),400);
    title('middle');

    subplot(4,1,4)
    auxplot(time,eeg.data,nsamp-400,400);
    title('end');

    in=input('Type "s" to save this tetrode, "n" to skip to next epoch:','s');
    oldepoch= -1;
    if length(in)==1
        switch lower(in)
        case 's'
            EEGselect{d(i)}{e(i)}= t(i);
        case 'n'
            oldepoch= e(i);
            oldday=d(i);
        end
    end
end
in=input('Press "y" to save EEG select to disk:','s');
if ~isempty(in) & in=='y'
    save EEGselect EEGselect
end

function auxplot(t,data,lo,steps)
hi=lo+steps-1;
plot(t(lo:hi),data(lo:hi));
xlabel('time [s]');
ylabel('theta amplitude [?]');
axis tight;

