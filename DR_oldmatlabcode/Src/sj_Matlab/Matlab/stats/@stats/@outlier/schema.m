function schema
%
% $Revision: 1.1.8.3 $	$Date: 2004/12/06 16:38:05 $
% Copyright 2003-2004 The MathWorks, Inc.

pk = findpackage('stats');

% Create a new class called outlier

c = schema.class(pk, 'outlier');

% Add properties
schema.prop(c, 'name', 'string');
p=schema.prop(c, 'dataset', 'string');

p=schema.prop(c, 'YLow', 'string');
p=schema.prop(c, 'YHigh', 'string');

% for these "equal" properties, "0" means "less/greater than or equal", 
% "1" means "less/greater than"
p=schema.prop(c, 'YLowLessEqual', 'double');
p=schema.prop(c, 'YHighGreaterEqual', 'double');

p=schema.prop(c, 'listeners', 'MATLAB array'); % place to store listeners
p.AccessFlags.Serialize = 'off';


