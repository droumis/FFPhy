% function [event] = extractevents(data, threshold, baseline, min_separation,
% 				   min_duration, allowstartend)
%	Extracts events where the data is above threshold for > min_duration 
%	samples.
%	
%	Adjacent events closer than min_separation samples are combined into a 
%	single event.
%	
%	Events are detected by threshold crossings.  The beginning of each
%	event is set to the closest sample preceeding the threshold crossing 
%	where the data went from below baseline to above baseline.  The end of
%	each event is set to the most distant sample following the threshold 
%	crossing where the signal is still above baseline.
%
%	allowstartend is 1 if events that begin or end above baseline should be
%	included and 0 if they should be excluded.
%
%	event is an 8xN matrix where each column represents an event.  The rows
%	are as follows:
%	row #		value
%	  1		start_index
%	  2		end_index
%	  3		peak_value
%	  4		peak_index
%	  5		total_area
%	  6		area_midpoint_index
%	  7		total_energy (e.g. sum squared values)
%	  8		energy_midpoint_index
%
%	Transposing the outputs of the function is recommended to maintaining
%	Frank Lab conventions, in which case the row #s above would refer to
%	column numbers.
