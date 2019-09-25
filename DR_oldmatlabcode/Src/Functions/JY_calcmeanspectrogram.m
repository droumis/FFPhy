function [ output_args ] = JY_calcmeanspectrogram(datadir,fileprefix,day,tet,reference)
%calculates the continuous spectrogram, epoch mean, and std spectra for each
%tetrode.
%   Based on Chronux package, parameters hard coded, check in code
%   datadir: directory where processed data is stored eg. /data18/jai/N2_/
%   fileprfix: name of rat
%   day: day to process
%   tetrodes: list of tetrodes to analyse eg. [1:5]
%   reference: provide a reference channel to add back referenced eeg
%   signal
%   Calculation is done for times during epochs specified in times.mat for
%   the day.
%   Derefencing takes account of time offset between tetrodes on different
%   DSPs.

dsz = '';
if (day < 10)
    dsz = '0';
end
%% load times.mat
filedir=datadir(1:end-2);
times=load(strcat(filedir,'/',fileprefix,dsz,num2str(day),'/times.mat'));
times=times.ranges(2:end,:);

d=day;
e=size(times,1);
t=max(tet);


%% load file
%load reference eeg tetrode
animaldir=strcat(datadir);

if ~isempty(reference)
    
    refeegstruct = loadeegstruct(animaldir, fileprefix,'eeg',day,[1:size(times,1)]',reference);
    
else
    refeegstruct=[];
end

% process epoch by epoch

d=day;
e=size(times,1);
t=max(tet);

continuousspectrograms=cell(d,e,t);
meandayspectra=cell(d,t);
stddayspectra=cell(d,t);

for epoch=1:size(times,1)
    
    eegstruct = loadeegstruct(animaldir, fileprefix,'eeg',day,epoch,tet);
    
    % go through each tetrode
    for ii=tet
        spectrograms={};     
        tetrodes=ii;
        
        % check if file exists
        tsz = '';
        if (tetrodes < 10)
            tsz = '0';
        end
        filename=strcat(datadir,'EEG/',fileprefix,'meaneegspectrograms',dsz,num2str(day),...
            '-',num2str(epoch),'-',tsz,num2str(tetrodes),'.mat');

        if exist(filename,'file')==0
            
            disp(sprintf('Epoch %s',num2str(epoch)));
            disp(sprintf('--Tetrode %s',num2str(tetrodes)));
            
            
            currtetstart=eegstruct{day}{epoch}{tetrodes}.starttime;
            currtetend=currtetstart+size(eegstruct{day}{epoch}{tetrodes}.data,1)/eegstruct{day}{epoch}{tetrodes}.samprate;
            spectrograms{d}{epoch}{tetrodes}.samprate=eegstruct{day}{epoch}{tetrodes}.samprate;
            % compare reference start and end times to find the times shared by
            % reference and current tetrode
            
            if ~isempty(reference)
                refeeg=refeegstruct{day}{epoch}{reference};
                refstart=refeeg.starttime;
                refend=refstart+size(refeeg.data,1)/refeeg.samprate;
                
                % if reference start time is earlier than current tetrode start time,
                % then start time is the later one
                if ~isequal(refstart, currtetstart)
                    if refstart<currtetstart
                        
                        adjrefstart=currtetstart;
                        adjrefstartind=ceil((currtetstart-refstart)/(1/refeeg.samprate));
                        adjrefsize=size(refeeg.data(adjrefstartind:end,1),1);
                        adjcurrstartind=1;
                        
                        % find the smaller of ref or current tetrode data length
                        newlength=min(adjrefsize,size(eegstruct{day}{epoch}{tetrodes}.data,1));
                    else
                        adjrefstart=refstart;
                        adjcurrtetstart=refstart;
                        adjrefstartind=1;
                        adjcurrstartind=ceil((refstart-currtetstart)/(1/refeeg.samprate));
                        adjcurrtetsize=size(eegstruct{day}{epoch}{tetrodes}.data(adjcurrstartind:end,1),1);
                        
                        % find the smaller of ref or current tetrode data length
                        newlength=min(adjcurrtetsize,size(refeeg.data,1));
                    end
                    
                    
                    adjrefdata=refeeg.data(adjrefstartind:newlength+adjrefstartind-1,1);
                    adjtetdata=eegstruct{day}{epoch}{tetrodes}.data(adjcurrstartind:newlength+adjcurrstartind-1,1);
                    
                    %add back reference
                    eegstruct{day}{epoch}{tetrodes}.data=adjrefdata+adjtetdata;
                    eegstruct{day}{epoch}{tetrodes}.starttime=adjrefstart;
                    
                    
                    spectrograms{d}{epoch}{tetrodes}.starttime=adjrefstart;
                    spectrograms{d}{epoch}{tetrodes}.data=adjrefdata+adjtetdata;
                    spectrograms{d}{epoch}{tetrodes}.unreferenced=eegstruct{day}{epoch}{tetrodes}.data;
                    
                else
                    
                    minlengthindex=min(size(eegstruct{day}{epoch}{tetrodes}.data,1),size(refeeg.data,1));
                    
                    spectrograms{d}{epoch}{tetrodes}.starttime=currtetstart;
                    spectrograms{d}{epoch}{tetrodes}.data=eegstruct{day}{epoch}{tetrodes}.data(1:minlengthindex)...
                        +refeeg.data(1:minlengthindex);
                    spectrograms{d}{epoch}{tetrodes}.samprate=eegstruct{day}{epoch}{tetrodes}.samprate;
                    spectrograms{d}{epoch}{tetrodes}.depth=eegstruct{day}{epoch}{tetrodes}.depth;
                    spectrograms{d}{epoch}{tetrodes}.unreferenced=eegstruct{day}{epoch}{tetrodes}.data;
                end
                
            else
                spectrograms{d}{epoch}{tetrodes}.starttime=currtetstart;
                spectrograms{d}{epoch}{tetrodes}.data=eegstruct{day}{epoch}{tetrodes}.data;
                spectrograms{d}{epoch}{tetrodes}.samprate=eegstruct{day}{epoch}{tetrodes}.samprate;
                spectrograms{d}{epoch}{tetrodes}.depth=eegstruct{day}{epoch}{tetrodes}.depth;
                spectrograms{d}{epoch}{tetrodes}.unreferenced=eegstruct{day}{epoch}{tetrodes}.data;
            end
            
            
            %% do spectral analysis
            
            % chronux params
%             movingwin = [125 10]/1000;
%             params.Fs = 1500;
%             params.err = [2 0.05];
%             params.fpass = [10 400];
%             params.tapers = [3 5];
%             
%             
%             dummy=[];
%             
%             [S_full,junkt,junkf,junkserr] = mtspecgramc(spectrograms{d}{epoch}{tetrodes}.data,movingwin,params);
%             dummy=[dummy;S_full];
%             continuousspectrograms{d}{epoch}{tetrodes}=S_full;
%             
%             meandayspectra{tetrodes}=mean(dummy,1);
%             stddayspectra{tetrodes}=std(dummy,1);
%             
%             % z-score the continuous spectrogram
%             
%             continuousspectrograms{d}{epoch}{tetrodes}=bsxfun(@minus,continuousspectrograms{d}{epoch}{tetrodes},meandayspectra{tetrodes});
%             continuousspectrograms{d}{epoch}{tetrodes}=bsxfun(@rdivide,continuousspectrograms{d}{epoch}{tetrodes},stddayspectra{tetrodes});
%             
%             spectrograms{d}{epoch}{tetrodes}.zscorecontinuous=continuousspectrograms{d}{epoch}{tetrodes};
%             spectrograms{d}{epoch}{tetrodes}.epochmean= meandayspectra{tetrodes};
%             spectrograms{d}{epoch}{tetrodes}.epochstd= stddayspectra{tetrodes};
            spectrograms{d}{epoch}{tetrodes}.reftet=reference;
            spectrograms{d}{epoch}{tetrodes}.descript=eegstruct{day}{epoch}{tetrodes}.descript;
            spectrograms{d}{epoch}{tetrodes}.fields=eegstruct{day}{epoch}{tetrodes}.fields;
            %save(filename,'-7.3');
            tsz = '';
            if (tetrodes < 10)
                tsz = '0';
            end

            meaneegspectrograms=spectrograms;

            
            filename=strcat(datadir,'EEG/',fileprefix,'meaneegspectrograms',dsz,num2str(day),...
                '-',num2str(epoch),'-',tsz,num2str(tetrodes),'.mat');
            
            save(filename,'meaneegspectrograms');
            disp(sprintf('---Saved %s', filename));
            
            
        else 
            fprintf('found %s \n',filename);
    
        end
    end
    
    
end


disp('.finished with computing spectrograms.')
end





