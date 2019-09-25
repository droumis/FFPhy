function out = parseFunctionCallString(fstring)
%fparams = parseFunctionCallString(fstring)
%
%Parses a string into parameters for a function call. 
%fstring sequence is '<function> FUNCNAME <argname> STRING <argval> VALUE <argname> ...'
%VALUE is evaluated, so numeric arrays, cell arrays, and strings are
%unserstood. For strings, use ''...''.
%fparams is a structure containing the function name and the call
%arguments.
%
%Example:
%fparams = parseFunctionCallString('<function> findpeople <argname> hair color <argval> ''brown'' <argname> age range <argval> [23:45]');
%

pat = '<[^>]+>'; %match pattern (anything that looks like this: <...>
controlArgs = regexp(fstring,pat,'match'); 
controlArgs = controlArgs(2:end);
functionArgs = regexp(fstring,pat,'split');
out.funcName = removeWhiteSpace(functionArgs{2},'all');
functionArgs = functionArgs(3:end);

if (length(functionArgs) ~= length(controlArgs))
    error('Incorrect syntax in function call string. Each <argname> must have a matching <argval>');
end
if mod(length(controlArgs),2)
    error('Incorrect syntax in function call string. Each <argname> must have a matching <argval>');
end
for i = 1:2:length(controlArgs)
    if ( (strcmp(lower(controlArgs{i+1}),'<argval>')) && strcmp(lower(controlArgs{i}),'<argname>'))                
        out.funcVarargin{i} = removeWhiteSpace(functionArgs{i},'ends');     
        out.funcVarargin{i+1} = eval(functionArgs{i+1});                           
    else
        error('Incorrect syntax in function call string. The call sequence is <function> FUNCNAME <argname> STRING <argval> VALUE <argname> ... ');
    end
end
%-------------------------------------------------------------
function strOut = removeWhiteSpace(strIn,mode);

if strcmp(mode,'all')
    strOut = regexprep(strIn,'[\s'']',''); %remove all whitespace
elseif strcmp(mode,'ends')
    strOut = regexprep(strIn,'^[\s]*',''); %removes whitespace only at ends 
    strOut = regexprep(strOut,'[\s]*$','');
end
    

