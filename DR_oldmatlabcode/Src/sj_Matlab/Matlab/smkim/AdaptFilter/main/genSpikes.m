function  spiketimes= genSpikes(data, model, generator)
%function  spiketimes= genSpikes(data, model, generator)
%
%     data:  timearray, posarray, phasearray
%     model: cpx, cpp, cpt, thetahat, thetainit
%
%	thetainit is the initial values of the control points (X first, then T)
%	KSStat is the KS statistics for the intervals
%	thetahat is the list of control point updates
%	spatialstats is a structure array with one element for each trajectory
%		containing the spatial statistics for that trajectory sampled
%		at the times given in statstimes
%	temporalstats is a structure containing the temporal statistics 
%		at the times given in statstimes
%
%	spiketimes is a 1xN array of the times for each spike [ms]
%	timearray is a 1xT list of times covering the entire period [ms]
%	posarray is a 1xT list of positions, one for each time 
%	epsx is a 1x4 array of the learning rates for the spatial spline
%	epst is a 1x4 array of the learning rates for the temporal spline
%	estinfo is a 1xT array of numbers indiciating which trajectory the
%		animal is on (start number at 0, -1 for points that should not be used in 
%		the estimation)
%	cpx is a cell array with P elements, one for each trajectory, where
%		each element lists the spatial control point locations for that 
%		trajectory
%	cpt is the list of temporal control point locations (in milliseconds)
%	niterations is the maximum number of iterations to run
%	statstimes is the list of times at which the mean, area, etc. for the
%		spatial and temporal spline are calculated. It should only be 
%		specified when there are five output arguments
%	statsegs is a Sx3 list containing, in each row, the trajectory number
%		the starting control point, and the ending control point for
%		each segment over which statistics should be calculated
