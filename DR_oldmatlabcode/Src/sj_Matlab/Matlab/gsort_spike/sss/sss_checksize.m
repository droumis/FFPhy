
%%% Check Size of raw spike data file in Continuous Recording Session, cut into 5 minute sections

fn='spike.dat'; Fs=32000; ch=[1:4];

for i=1:15   
    t1=300*(i-1); t2=300*i;    
    data=labview_loadspike4(fn,4,Fs,t1,t2);  %% NEED TO READ ALL 8 channels
    data=data(ch,:); data=data';
    
    %file=['tempdatan' num2str(i)];
    %save(file, 'data');
    if size(data,1)<Fs*300
        i
        break
    end
end