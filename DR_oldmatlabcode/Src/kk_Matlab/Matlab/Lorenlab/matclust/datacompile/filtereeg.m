function filteredeegstruct = filtereeg(raweegstruct, filterstruct)
%
%   filtereeg(raweegstruct, filterstruct)
%
% Applies the specified zero-phase finite impulse response FILTER to the 
% the specified RAWEEG structure. Returns a FILTEREDEEG structure whose 
% data field contains the filtered amplitude, instantaneous phase 
% (extracted by Hilbert transform), and envelope magnitude (extracted by 
% Hilbert transform). 
%
% Note that the sampling rate of the FILTER must match the sampling rate 
% of the eeg data. Also, note that FILTER must have a finite impulse 
% response which is symmetric about zero and has an odd number of 
% coefficients. FILTER is a special structure with fields 'descript', 
% 'kernel' and 'samprate'.
%
% Finally, be aware that the very beginning and end time segments of the
% filtered output are not to be trusted because the required samples are
% missing from the input.
%
% written by smk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERROR CHECKING


% Check that the supplied raweeg data is a valid column vector
if isempty(raweegstruct.data)
    error('eeg data is empty')
elseif length(raweegstruct.data) ~= numel(raweegstruct.data)
    error('eeg data must be a column vector')
else
    [numrows, numcols] = size(raweegstruct.data);
    if numrows==1
        error('eeg data must be a column vector')
    end
end

% Check that the filter kernel is a non-empty vector
if isempty(filterstruct.kernel)
    error('filter kernel is empty')
elseif length(filterstruct.kernel) ~= numel(filterstruct.kernel)
    error('filter kernel must be a vector')
end

% Error checking on the samprates
samprate_ratio = raweegstruct.samprate/filterstruct.samprate;

% Check that sampling rates of the filter and the data match to within some
% very small epsilon; warn the user if this is not the case
if abs(samprate_ratio - 1) > 1e-3
    warning('filterstruct.samprate does not match eegstruct.samprate');
end

% Check that filterstruct.samprate does not exceed eegstruct.samprate
if (samprate_ratio - 1) < -1e-3
    disp(samprate_ratio);
    error('filterstruct.samprate must not exceed eegstruct.samprate');
end

% Check that samprate_ratio is close to an integer
if abs(round(samprate_ratio) - samprate_ratio) > 1e-3
    error('filterstruct.samprate must be very close to an integer fraction of eegstruct.samprate');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW FOR THE ACTUAL BUSINESS

% copy the entire structure
filteredeegstruct = raweegstruct;

% !!! NOTE that filtfilt does not handle NaN entries. If the data contain NaN
% values, we need to split the data into continuous segments of non-NaN values,
% apply filtfilt to each separately, and then join then back together with 
% the correct padding intervals of NaN values
nanidxs = find(isnan(raweegstruct.data))';
nanidxs = [ 0 nanidxs (length(raweegstruct.data)+1) ];
nanintervals = find(diff(nanidxs)>1);

% allocate array to hold the filtered output, initialized with NaN
filteredeegstruct.data = NaN(length(raweegstruct.data),3);

% filter each continuous interval of non-NaN values and copy into the 
% appropriate section of the output array (leaving NaN values in between)
for i = 1:length(nanintervals)
    idx_start = nanidxs(nanintervals(i)) + 1;
    idx_end = nanidxs(nanintervals(i) + 1) - 1;
    % now filter each segment separately
    if round(samprate_ratio) == 1
        % if filterstruct.samprate exactly matches eegstruct.samprate, we 
        % filter over the dataset user the Signal Processing toolbox 
        % filtfilt function
        try
            filteredeegstruct.data(idx_start:idx_end,1) = filtfilt( ...
                filterstruct.kernel,[1],...
                raweegstruct.data(idx_start:idx_end)); 
            % apply Hilbert transform
            temp = hilbert(filteredeegstruct.data(idx_start:idx_end,1));
            filteredeegstruct.data(idx_start:idx_end,2) = angle(temp);
            filteredeegstruct.data(idx_start:idx_end,3) = abs(temp);
        catch
            % the data segment may be too short to filter over,
            % in which case we skip it
        end
    else
        % otherwise, we need to decimate, filter, and then interpolate
        % Two things to note: First, it is very important to specify 'fir' 
        % (finite impulse resopnse filter) in decimate; the default 
        % Chebyshev filter introduces a phase delay (very naughty!). 
        % Second, decimate yields a signal whose number of samples is 
        % ceil(length(raweegstruct.data)/samprate_ratio); therefore, 
        % unless mod(length(raweegstruct.data),round(samprate_ratio))==0, 
        % the subsequent interp will yield an output of slightly different 
        % length than the original data vector. If I did this all correctly, 
        % the starttime should coincide and the exra samples at the end can 
        % be safely ignored (I hope)
        try
            temp = decimate( ...
                raweegstruct.data(idx_start:idx_end),round(samprate_ratio), ...
                400,'fir');
            temp = filtfilt(filterstruct.kernel,[1],temp);
            temp = interp(temp,round(samprate_ratio));
            % interp may give us a few extra points at the end, which we 
            % need to trim off
            temp = temp(1:(idx_end-idx_start+1));
            filteredeegstruct.data(idx_start:idx_end,1) = temp;
            % apply Hilbert transform
            temp = hilbert(filteredeegstruct.data(idx_start:idx_end,1));
            filteredeegstruct.data(idx_start:idx_end,2) = angle(temp);
            filteredeegstruct.data(idx_start:idx_end,3) = abs(temp);
        catch
            % the data segment may be too short to filter over,
            % in which case we skip it
        end
    end
end


% overwrite the fields field to reflect the contents of the three columns
filteredeegstruct.fields = ...
'filtered_amplitude1 instantaneous_phase2 envelope_magnitude3';

% add in a copy of the filter kernel
filteredeegstruct.filterkernel = filterstruct.kernel;
filteredeegstruct.filtersamprate = filterstruct.samprate;

% finally edit the descript field of this cell to reflect that this is 
% *filtered* data, not the original raw eeg trace
if round(samprate_ratio)==1
    filteredeegstruct.descript = ...
    strvcat(raweegstruct.descript,'filtered through: ', ...
    [blanks(size(filterstruct.descript,1))' filterstruct.descript]);
else
    filteredeegstruct.descript = ...
    strvcat(raweegstruct.descript,'decimated, filtered through: ', ...
    [blanks(size(filterstruct.descript,1))' filterstruct.descript], ...
    'and interpolated back to original samprate');
end
