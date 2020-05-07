function [eegvar] = iseegvar(varname)
% eegvar = iseegvar(varname)
% returns 1 if varname is the name of an eeg variable and 0 otherwise
% 
% eeg variables are
% 	'eeg'
%	'theta'
%	'gamma'
%	'ripple'
switch (varname)
case {'eeg', 'theta', 'gamma', 'ripple','eeggnd','thetagnd','gammagnd','ripplegnd', 'slowgamma', 'fastgamma', 'lowgamma', 'highgamma'}
    eegvar = 1;
otherwise
    eegvar = 0;
end

