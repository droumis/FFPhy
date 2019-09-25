function [data] = geteegdata(eegstruct, datatype)
% function [data] = geteegdata(eegstruct, datatype)
%     Given data type, returns the value of the eeg data
%	datatype is one of the following:
%	'amp'	- returns the amplitude of the signal
%	'phase'	- returns the phase of the signal, adjusting for the 10000x
%		   phase multiplier
%	'env'	- returns the magnatude of the envelope
switch (datatype)
case 'amp'
    data = double(eegstruct.data(:,1));
case 'phase'
    data = double(eegstruct.data(:,2)) / 10000;
case 'env'
    data = double(eegstruct.data(:,3));
otherwise
    disp('geteegdata: unknown datatype');
end

