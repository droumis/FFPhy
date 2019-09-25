function [out] = structure_group_combine(f, field_1, field_2, varargin)
% function [out] = structure_group_combine(f, field_1,field_2, varargin)
% A structure version of numericgroupcombine.  Collapses filter fields across animals.  
% field_1  Is the string 'field_name' referring to the output field.  
%          ie, if field_1 is set to 'output' this function will combine contents, across animals, for
%          f(an).output. 
% ***** if field_1 = 'epochs', a particular for loop is invoked and field_2 will not be used ***
% 
% field_2  String with the 'field_name_2' of the subfield.  
%          It can be empty '' if not needed, that is, if f(an).(field_1) is numeric.  
%          Otherwise it is the 'field_name_2' of a field in an output structure.
%          ie, if field_1 is assigned as above and field_2 is set to 'peak_rates'
%          you will be combining contents, across animals, for
%          f(an).output{epochID}.peak_rates
%          WARNING: I HAVEN'T TRIED THIS YET WITH EMPTY FIELD_2
% NOTE: expects that the data to be combined are numeric.  
% 
%
% OPTIONS
% emptyfield  [empty field contents]  Send in whatever you'd like stored if it's empty.
%                                     Default empty vector. This is mostly useless, it's
%                                     only employed if there are no good data in the entire filter.
%
% anathe JUL2009

empty_contents = [];
stacked_fields = [];

for option = 1:length(varargin)-1
    switch lower(varargin{option})
    case 'emptyfield'
	empty_contents = varargin{option+1}
    otherwise
	error (['Varargin ', varargin{option}, ' not defined.'])
    end
end

if isequal(field_1,'epochs')
    for an = 1:length(f)
	if ~isempty(f(an))

	    for epoch_filter = 1:length(f(an).(field_1))
		if ~isempty(f(an).(field_1){epoch_filter})
		    stacked_fields = [stacked_fields; f(an).(field_1){epoch_filter}];
		end
	    end

	end
    end
    out = stacked_fields;
    return
end


for an = 1:length(f)
    if ~isempty(f(an))
	
	for epoch_filter = 1:length(f(an).(field_1))
	    if ~isempty(f(an).(field_1){epoch_filter})

		for day_epoch_pair = 1:length(f(an).(field_1){epoch_filter})
		    if ~isempty(f(an).(field_1){epoch_filter})
			
			if  isempty(field_2)
			    field_contents = f(an).(field_1){epoch_filter}(day_epoch_pair);
			else
			    field_contents = f(an).(field_1){epoch_filter}(day_epoch_pair).(field_2);
			end
		    end
		    if isempty(find(~isnan(field_contents)))
			continue;
		    end
		    stacked_fields = [stacked_fields; field_contents];
		    clear field_contents;
		end

	    end
	end

    end
end

if isempty(stacked_fields)
    stacked_fields = empty_contents
end

out = stacked_fields;	   


