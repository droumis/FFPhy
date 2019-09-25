function phase = estimate_oscillatory_phase(continuous,method,expected_period,varargin)
%ESTIMATE_OSCILLATORY_PHASE Estimate phase of oscillations in a bandpass-filtered signal
%
%   This function provides methods to estimate the phase of an oscillation
%   without relying on an analytic representation of the signal as a narrowband
%   sinusoid. These methods were applied to theta oscillations, which have a
%   pronounced sawtooth-like waveform assymmetry (Siapas et al., 2005). In
%   principle, the same methods could also be applied to other oscillatory
%   signals which can not be modeled as narrowband modulated sinusoids, such as
%   neocortical spindle oscillations or even periodically-driven field
%   responses.
%   
%   PHASE = ESTIMATE_OSCILLATORY_PHASE(CONTINUOUS,METHOD,EXPECTED_PERIOD) takes
%   a CONTINUOUS data struct and returns a PHASE data struct which inherits all
%   of the same fields as the original, except that the 'samples' field contains
%   complex values whose angle is the estimated phase and whose magnitude is 1.
%   If the 'samples' field of the input CONTINUOUS is complex, only its real
%   part is considered. METHOD can be one of the following strings:
%     'minima': Minima in the waveform are assigned phase of +pi == -pi, and
%         between minima the phase is linearly interpolated.
%     'maxima': Maxima in the waveform are assigned phase of 0, and between
%         maxima the phase is linearly interpolated.
%     'extrema': Maxima are assigned phase of +pi and minima are assigned phase
%         of 0, and over each half-cycle the phase is linearly interpolated
%     'up': zero-crossings from negative to positive are assigned phase of
%         -pi/2, and between these upwards crossings the phase is linearly
%         interpolated
%     'down': Zero-crossings from positive to negative are assigned phase of
%         +pi/2, and between these downwards crossings the phase is linearly
%         interpolated
%     'zero_crossings': Zero-crossings in both directions are assigned the
%         respective phase -pi/2 and +pi/2, and over each half-cycle the phase
%         is linearly interpolated
%   METHOD is appended to the PHASE output struct as a
%   'phase_estimation_method' field. EXPECTED_PERIOD is a finite positive real
%   floating-point scalar which specifies the approximate mean period (in
%   seconds) of the oscillation. This parameter can be estimated from the
%   spectrum of the LFP or by visual inspection of the trace.
%
%   For the phase estimation to be effective, the input signal in CONTINUOUS
%   should be pre-processed to remove noise and DC offset. (For theta
%   oscillations, the recommended pre-processing is to smooth with a 4.5-40 Hz
%   bandpass filter. The filter should be applied in a way that does not
%   introduce phase distortion.) Alternatively, instead of pre-processing the
%   input signal, you can pass a CONTINUOUS data struct containing raw signal
%   along with a filter specification object (and optionally, a filter object
%   too):
%
%   ESTIMATE_OSCILLATORY_PHASE(CONTINUOUS,METHOD,D)
%   ESTIMATE_OSCILLATORY_PHASE(CONTINUOUS,METHOD,D,HD)
%
%   See the help for FILTER_CONTINUOUS (written by SMK) for further information
%   on these filtering arguments D and HD.
%
%   Reference: 
%
%   Siapas A., Lubenov E., Wilson M. (2005) Prefrontal phase locking to
%   hippocampal theta oscillations. _Neuron_ 46:141-151.
%
%Depends on:
%   IS_CONTINUOUS (written by SMK)
%   INTERP1_ANGLE (written by SMK)
%   FILTER_CONTINUOUS (written by SMK)
%
%Written by SMK, 2009 November 20.
%

  TS_PER_SEC = 1e4;

  if (exist('is_continuous') ~= 2)
    error(['ESTIMATE_OSCILLATORY_PHASE depends on m-file IS_CONTINUOUS ' ...
        '(written by smk)']);
  end
  if (exist('interp1_angle') ~= 2)
    error(['ESTIMATE_OSCILLATORY_PHASE depends on m-file INTERP1_ANGLE ' ...
        '(written by SMK)']);
  end
  if (exist('filter_continuous') ~= 2)
    error(['ESTIMATE_OSCILLATORY_PHASE depends on m-file FILTER_CONTINUOUS ' ...
        '(written by SMK)']);
  end

  if ~is_continuous(continuous)
    error(['CONTINUOUS does not appear to be a valid continuous data ' ...
        'struct array']);
  end
  if ~ischar(method)
    error('METHOD must be one of the recognized method names');
  end
  if ~any(strcmp(method, ...
      {'minima','maxima','extrema','up','down','zero_crossings'}))
    error('%s is not a recognized estimation method name',method);
  end

  if ~iscalar(expected_period) || ~isreal(expected_period) || ...
      ~isfloat(expected_period) || ~isfinite(expected_period) || ...
      ~(expected_period > 0)
    error(['EXPECTED_PERIOD must be a positive finite real floating-point' ...
        ' scalar']);
  end

  if (length(varargin) > 1)
    error('Too many additional arguments');
  elseif (length(varagin) > 0)
    % Replace original signal with filtered signal
    try
      continuous = filter_continuous(continuous,varargin{:});
    catch
      error(['FILTER_CONTINUOUS failed. Check that filter argument is a ' ...
          'MATLAB Signal Processing Toolbox ''dfilt'' class object.']);
    end
  end

  % Inherit meta-data fields
  phase = rmfield(continuous,'samples');
  % Add the method as a field
  [phase(:).phase_estimation_method] = deal(method);
  for i = 1:numel(phase)
    % Convert uint32 timestamps to double seconds
    ti = double(phase(i).timestamp)/TS_PER_SEC;
    try
      [maxima, down, minima, up] = ...
          find_cycles(ti,double(real(continuous(i).samples)),expected_period);
    catch
      error(['Error in the find_cycles routine with expected period of %f ' ...
          'seconds'],expected_period);
    end
    switch method
    case 'minima'
      t = minima;
      a = 2*pi*(1:numel(minima))' - pi;
    case 'maxima'
      t = maxima;
      a = 2*pi*(1:numel(maxima))'
    case 'extrema'
      t = [ minima; maxima ];
      if (minima(1) < maxima(1))
        a = [ 2*pi*(1:numel(minima))' - pi; 2*pi*(1:numel(maxima))' ]; 
      else
        a = [ 2*pi*(1:numel(minima))' + pi; 2*pi*(1:numel(maxima))' ]; 
      end
      [t, idx] = sort(t,1,'ascend');
      a = a(idx);
    case 'up'
      t = up;
      a = 2*pi*(1:numel(up))' - pi/2;
    case 'down'
      t = down;
      a = 2*pi*(1:numel(down))' + pi/2;
    case 'zero_crossings'
      t = [ up; down ];
      if (up(1) < down(1))
        a = [ 2*pi*(1:numel(up))' - pi/2; 2*pi*(1:numel(down))' + pi/2 ];
      else
        a = [ 2*pi*(1:numel(up))' + 3*pi/2; 2*pi*(1:numel(down))' + pi/2 ];
      end
      [t, idx] = sort(t,1,'ascend');
      a = a(idx);
    otherwise
      error(['Method %s is unrecognized; this error should have been ' ...
          'caught earlier during input checking.'],method);
    end
    % cos and sin accumulate errors when their arguments are large, so we first
    % wrap mod (2*pi)
    ai = recenter_angle(interp1(t,a,t_i,'linear'),0);
    phase(i).samples = single(complex(cos(ai),sin(ai)));
  end

end % end main function ESTIMATE_OSCILLATORY_PHASE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given time series (t,x), find complete oscillatory cycles and the occurrence
% times of the features [maxima, down, minima, up]. The four outputs are column
% vectors that may not all the same size.
function [maxima, down, minima, up] = find_cycles(t,x,period)


end % end function find_cycles


% Find the times of zero-crossings, given two bounding points at each location.
% This function takes vector inputs.
function ti = crossings_interp(t1,t2,x1,x2)
  assert(all(t2 > t1) && all(sign(x1) ~= sign(x2)));

  ti = (t1 .* x2 - t2 .* x1) ./ (x2 - x1);
end


% Find the times of maxima/minima, given three bounding points at each location.
% This function takes vector inputs.
function ti = extrema_interp(t1,t2,t3,x1,x2,x3)
  assert(all(t2 > t1) && all(t3 > t2) && all(sign(x2 - x1) ~= sign(x3 - x2)));

  a = t3 .* (x2 - x1) + t2 .* (x1 - x3) + t1 .* (x3 - x2);
  b = t3.^2 .* (x1 - x2) + t2.^2 .* (x3 - x1) + t1.^2 .* (x2 - x3);
  ti = -b ./ (2*a);
end


