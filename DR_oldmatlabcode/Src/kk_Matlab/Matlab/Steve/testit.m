function position = estimate_position(rawpos,varargin)
%ESTIMATE_POSITION Estimate head direction, smoothed position and velocity from raw video pixel coordinates of head-tracking LEDs.
%
%   POSITION = ESTIMATE_POSITION(RAWPOS,SMOOTHING_PARAMS) estimates
%   head direction, smoothed position and velocity from the coordinates of the
%   head-tracking LEDs in the RAWPOS.
%
%   Next, the line segment is drawn connecting the positions of front and back
%   LEDs in the real-world coordinate frame (e.g. units of centimeters). The
%   direction of this line segment is filtered to estimate head direction. A
%   point along this line segment (between the two LEDs) is chosen and its
%   location is smoothed with LOESS regression. Stopping times are identified
%   by applying a repeated running median smoother to the smoothed position to
%   identify stops (note that "stopped" does not necessarily mean immobile or
%   asleep; whisking, grooming, etc. may be occuring).
%
%   Optional arguments:
%
%   centimeters_per_pixel
%     Parameter for scaling pixels. Default is 1 cm / pixel.
%
%   front_back_marker_weights (applicable only when both back and front markers
%   are tracked)
%     Parameter for the relative weighting of the front LEDs versus back LEDs
%     when estimating "true" head position. This is a two-element vector of
%     positive reals. Typical value is [1 0] (i.e. rely on the front LEDs).
%
%   loess_halfwidth*
%     The halfwidth (expressed in seconds) of the sliding window on which
%     LOESS is done. This is the most important parameter for estimating
%     smooth position and velocity. Typical value is 0.75 (750 ms).
%
%   num_loess_iterations*
%     Number of reweighting iterations to apply LOESS for robust rejection of
%     outliers. Unless the data are especially noisy, the estimate practically
%     converges after 3 or 4 iterations.
%
%   rrm_halfwidths*
%     Vector of halfwidths (in seconds) of sliding windows in which running
%     medians are successively computed for the repeated running median
%     smoother. Data are smoothed with a running median of halfwidth
%     rrm_halfwidths(1), and then this result is smoothed with a running median
%     of rrm_halfwidths(2), etc. Window sizes must be non-increasing:
%     rrm_halfwidths(1) >= rrm_halfwidths(2) >= ... >= rrm_halfwidths(end).
%     Typical value is [0.1 0.1 0.1].
%
%   epsilon*
%     The "closeness" threshold (in cm) for determining when the subject is
%     immobile. This should be a very small value, as small as or smaller than
%     the pixel resolution. Typical value is 0.5.
%
%   min_stop_duration*
%     Define the minimum duration (in seconds) of consecutive "close" position
%     samples that must occur for an interval to qualify as a stop. Typical
%     value is 0.5 (500 ms).
%
%   min_front_back_separation (applicable only when both back and front markers
%   are tracked)
%     When the pitch of the rat's head is excessive, the front and back lights
%     are not distinguishable any more; therefore, head direction estimates are
%     unreliable at these times. This (float) parameter is the minimum apparent
%     front-back separation (in cm) required for estimating head direction.
%     Obviously, the value depends on the size of the LED tracker boom.
%     Default value is 2 cm.
%
%   *The smoothing parameters that are marked with asterisk are described in
%   Hen I., Sakov A., Kafkafi N., Golani I., Benjamini Y. (2004) The dynamics
%   of spatial behavior: how can robust smoothing techniques help? _Journal of
%   Neuroscience Methods_ 133: 161-172.
%
%   For explanation of regression model that is used to estimate head
%   direction, see section 6.4 of Fisher N.I. (1996)  _Statistical Analysis of
%   Circular Data_. Cambridge UP.
%
%Depends on:
%   TFORMFWD (MATLAB Image Processing Toolbox)
%   TFORMWD (MATLAB Image Processing Toolbox)
%   PDIST (MATLAB Statistics Toolbox)
%   FIND_NEARBY_MEX (written by smk)
%   BISQUARE (written by smk)
%   TRICUBE (writtten by smk)
%
%Written by smk (2009 November 1)
%
centimeters_per_pixel = 1;
front_back_marker_weights = [1 0];
loess_halfwidth = 0.75;
num_loess_iterations = 4;
rrm_halfwidths = [0.1 0.1 0.1];
epsilon = 0.5;
min_stop_duration = 0.5;
min_front_back_separation = 2;
TS_PER_SEC = 1e4;

[otherArgs,params] = procOptions(varargin);

keyboard
