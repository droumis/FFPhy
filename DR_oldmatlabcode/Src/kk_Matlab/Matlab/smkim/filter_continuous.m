function filtered = filter_continuous(continuous,Hd,varargin)
%FILTER_CONTINUOUS Filter continuous data.
%
%   FILTERED = FILTER_CONTINUOUS(CONTINUOUS,HD) filters the signal in the
%   continuous data struct CONTINUOUS and returns a matching FILTERED struct.
%   Samples in CONTINUOUS must have evenly-spaced timestamps; an error will be
%   raised if there are departures from near-constant sampling rate. FILTERED
%   exactly inherits all fields of CONTINUOUS, except that its 'samples' field
%   contains a filtered signal.
%
%   HD must be a 'dfilt' class object which is applied to the signal. FILTERED
%   contains a copy of this object as its 'filter_object' field. It is
%   the caller's responsibility to make sure that Hd is correctly specified for
%   the sampling rate of the original continuous data.
%
%   FILTER_CONTINUOUS(CONTINUOUS,HD,'phase') also compute instantaneous phase
%   using the Hilbert transform. The result is meaningful only when using a
%   bandpass filter that matches a narrowband oscillation in the signal.
%
%   The return struct FILTERED has the following fields:
%     subject: inherited from CONTINUOUS
%     day: inherited from CONTINUOUS
%     epoch: inherited from CONTINUOUS
%     electrode: inherited from CONTINUOUS
%     channel: inherited from CONTINUOUS
%     depth: inherited from CONTINUOUS
%     region: inherited from CONTINUOUS
%     reference: inherited from CONTINUOUS
%     Fs: the nominal sampling rate (double scalar)
%     timestamp: monotonically-increasing vector of uint32 timestamps
%     samples: vector of single values produced by the filter
%     (phase: vector of double values, equal to instantaneous phase expressed in
%         radians on the interval [-pi,+pi]. this field is included only when
%         the user calls the function with the 'phase' option)
%     source: inherited from CONTINUOUS
%     filter_object: a copy of the dfilt object HD
%
%Depends on:
%   IS_CONTINUOUS (written by smk)
%   HILBERT (MATLAB Signal Processing toolbox)
%   DFILT.* and FDESIGN.* package classes (MATLAB Signal Processing toolbox)
%   FILTFILTHD (written by Malcolm Lidierth, MATLAB Central File ID #17061; see 
%       http://www.mathworks.com/matlabcentral/fileexchange/17061-filtfilthd)
%
%Written by SMK, 2009 February 05.
%

if (exist('is_continuous') ~= 2)
  error('FILTER_CONTINUOUS depends on m-file IS_CONTINUOUS (written by smk)');
end
if (exist('hilbert') ~= 2)
  error(['FILTER_CONTINUOUS depends on m-file HILBERT ' ...
      '(MATLAB Signal Processing toolbox)']);
end
if (exist('filtfilthd') ~= 2)
  error(['FILTER_CONTINUOUS depends on m-file FILTFILTHD ' ...
      '(written by Malcolm Lidierth, MATLAB Central File ID #17061)']);
end

if ~is_continuous(continuous)
  error(['CONTINUOUS does not appear to be a valid continuous data ' ...
      'struct array']);
end
% TODO: Use meta-classes to interrogate class identity
if ~strncmp('dfilt.',class(Hd),6)
  error('HD must be a MATLAB Signal Processing Toolbox dfilt object');
end

if isempty(varargin)
  compute_phase = false;
  compute_envelope = false;  
elseif ~iscellstr(varargin) || ...
    ~all(strcmp('phase',varargin) | strcmp('envelope',varargin))
  error('Optional arguments must be the string ''phase'' and/or ''envelope''');
else
  compute_phase = any(strcmp('phase',varargin));
  compute_envelope = any(strcmp('envelope',varargin));
end

% Inherit meta-data fields
filtered = rmfield(continuous,'samples');
% Add the filter specification as a field
[filtered(:).filter_object] = deal(Hd);
for i = 1:numel(continuous)
  % Perform the computations on double precision for numerical stability
  try
    signal = filtfilthd(Hd,double(continuous(i).samples),'reflect');
    assert(all(isfinite(signal(:))));
  catch
    error(['Could not filter, or filter application resulted in ' ...
        'numerical instability']);
  end
  filtered(i).samples = cast(signal,class(continuous(i).samples));
  if (compute_phase || compute_envelope)
    analytic_signal = hilbert(signal);
    if compute_phase
      filtered(i).phase = angle(analytic_signal);
    end
    if compute_envelope
      filtered(i).envelope = abs(analytic_signal);
    end
  end
end


