function [zspectrogram,times,frequencies] = waveletterz(data,Fs,scale,mother,varargin)

% z-scores waveletter() spectrogram to all of data.

% Accepts either a vector of eeg data OR an eegstruct w/ .data fields
    % if an eegstruct, include DAY, EPOCH, and TETRODE in last 3 args
    
    % e.g. [S,t,f]=waveletterz(eegstruct,Fs,scale,mother,d,e,t)
    
% Note that if data is eegstruct w/ multiple tetrodes, then z-scores respective tetrode
% Note that longer epochs will have more lower frequencies represented --
    % these extra lower frequencies are truncated in order to perform the
    % z-scoring across epochs.
    
days = varargin{1};      
epochs = varargin{2};
tetrodes = varargin{3};

if isnumeric(data)
    [spectrogram,times,frequencies,~] = waveletter(data,Fs,scale,mother);
    meanspectrum=mean(spectrogram,2);
    stdspectrum=std(spectrogram,0,2);
    zspectrogram=bsxfun(@minus,spectrogram,meanspectrum);
    zspectrogram=bsxfun(@rdivide,zspectrogram,stdspectrum);
elseif iscell(data)
    for t=tetrodes
        dummy{t}=[];
        minepochlength=inf;
        frequencies=[];
        for d=days
            for e=epochs
                if ~isempty(data{d}{e}{t}.data)
                    [zspectrogram{d}{e}{t},~,f,~]=waveletter(data{d}{e}{t}.data,Fs,scale,mother);
                    if length(data{d}{e}{t}.data) < minepochlength
                        frequencies=f;
                        minepochlength=length(data{d}{e}{t}.data);
                    end
                    dummy{t}=[dummy{t} zspectrogram{d}{e}{t}(1:length(frequencies),:)];
                end
            end
        end
    end
    if ~isempty(dummy)                  % z-scoring with lumped epochs
        for t=tetrodes
            if ~isempty(dummy{t})
                meanspectrum{t}=mean(dummy{t},2);
                stdspectrum{t}=std(dummy{t},0,2);
            for d=days
                for e=epochs
                        zspectrogram{d}{e}{t}=zspectrogram{d}{e}{t}(1:length(frequencies),:);
                        zspectrogram{d}{e}{t}=bsxfun(@minus,zspectrogram{d}{e}{t},meanspectrum{t});
                        zspectrogram{d}{e}{t}=bsxfun(@rdivide,zspectrogram{d}{e}{t},stdspectrum{t});

                end
            end
            end
        end
        times=NaN;
    else
        error('no data found!')
    end
else
    error('data must be either an eegstruct{d}{e}{t} with .data field or raw eeg');
end



end

