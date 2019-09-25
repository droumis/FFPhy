function [zspectrogram,times,frequencies] = mtspecgramcz(data,movingwin,params,varargin)

% z-scores mtspecgram() spectrogram to all inputted data

% Accepts either a vector of eeg data OR an eegstruct w/ .data fields
    % if an eegstruct, include DAY, EPOCH, and TETRODE in last 3 args

    % e.g. [S,t,f]=waveletterz(eegstruct,Fs,scale,mother,d,e,t)
    
% Note that if data is eegstruct w/ multiple tetrodes, then z-scores
% respective tetrode.

% Also note that if data is an eegstruct, times is outputted for each epoch.

days = varargin{1};      
epochs = varargin{2};
tetrodes = varargin{3};

if isnumeric(data)
    [spectrogram,times,frequencies,~] = mtspecgramc(data,movingwin,params);
    spectrogram=spectrogram';
    meanspectrum=mean(spectrogram,2);
    stdspectrum=std(spectrogram,0,2);
    zspectrogram=bsxfun(@minus,spectrogram,meanspectrum);
    zspectrogram=bsxfun(@rdivide,zspectrogram,stdspectrum);
elseif iscell(data)
    times=cell(size(data));         % since movingwin dictates times, useful to store times in its own vector
    dummy=cell(1,tetrodes(end));
    for t=tetrodes
        for d=days
            for e=epochs
                if ~isempty(data{d}{e}{t}.data)
                    [zspectrogram{d}{e}{t},times{d}{e}{t},frequencies]=mtspecgramc(data{d}{e}{t}.data,movingwin,params);
                    dummy{t}=[dummy{t} zspectrogram{d}{e}{t}'];
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
                        zspectrogram{d}{e}{t}=bsxfun(@minus,zspectrogram{d}{e}{t},meanspectrum{t}');
                        zspectrogram{d}{e}{t}=bsxfun(@rdivide,zspectrogram{d}{e}{t},stdspectrum{t}');
                end
            end
            end
        end
    else
        error('no data found!')
    end
else
    error('data must be either an eegstruct{d}{e}{t} with .data field or raw eeg');
end



end

