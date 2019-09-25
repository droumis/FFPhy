function[f] = add_filter_info(f, field_name_string, field_contents)
% function[f] = add_filter_animal_info(f, field_name_string, field_contents)
% 
% Tacks a user-defined field onto the filter, for each animal.
% I'm using it to store the function call information for preliminary analyses,
% ie if I run a peakrate function using the filter, and use the result as a filter input
% to subsequent filter functions, I'd store the peakrate function call as the 
% contents of a peakrate_function field.  
%
% anathe may2009

for i = 1:length(f)
    f(i).(field_name_string) = field_contents;
end
    


