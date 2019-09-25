function r = rules (label, start, finish, cues, conditions, blocks, responses, trials)

% @rules/rules.m
% creates an object of class @rules; this object is passed to SPIKEMATRIX
% to select time windows from a @cortex object.
%
% Syntax 1:
% ROBJ = RULES (LABEL, START, FINISH, CUES, CONDITIONS, BLOCKS, RESPONSES, TRIALS)
% creates a one-rule object and assigns the values passed as arguments
% to the corresponding fields.
%
% LABEL:         the name of the analysis;
%                defaults to 'noname' if an empty string is passed;
% START, FINISH: the offset in ms from the cue event;
% CUES:          an array of encodes to time-lock the analysis on;
% CONDITIONS:    see below
% BLOCKS:        see below
% RESPONSES:     see below
% TRIALS:        arrays of values to consider for the analysis; use -1 to include
%                all (e.g., if r.blocks == -1 then all blocks are included in the 
%                analysis performed by GETSPIKEARRAY).
%
% Syntax 2:
% ROBJ = RULES
% creates a one-rule object and sets its fields to the following default values:
% LABEL = 'noname', START = -50, FINISH = 250, CUES = 23, RESPONSES = 0;
% CONDITIONS = BLOCKS = TRIALS = -1.
%
% Syntax 3:
% ROBJ = RULES (GROUPSFILE)
% creates a multiple rules/groups object based on
% an ASCII file of name GROUPSFILE
%
% Once the object is created, all the fields
% listed above can be read/modified, e.g.:
%
% » robj = rules;
% » robj.conditions = [4:8,12:16];
% » current_analysis = robj.label;
%
% Use the property SIZE to get the number of rules and groups
% contained in an object:
%
% » robj = rules ('groups.txt');
% » s = robj.size;
% » last_group = length (s);
% » last_rule  = s (last_group);
% » disp (robj (last_group, last_rule).conditions);
%
% See also: CORTEX, SPIKEMATRIX
% Last modified: 26 Nov 98

% if called with no arguments, a single 'rules' object is created
% with all default values

switch nargin
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0  % default object
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   r.st.label  = 'noname';
   r.st.start  = -50;
   r.st.finish = 250;
   r.st.cues   =  23;
   r.st.conds  =  -1;
   r.st.blocks =  -1;
   r.st.resps  =   0;
   r.st.trials =  -1;
   r.size      =   1;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 8  % then use them all
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if isempty (label)
      r.st.label = 'noname';
   else
      r.st.label = label;
   end
   if isempty (start)
      r.st.start = -50;
   else
      r.st.start = start;
   end
   if isempty (finish)
      r.st.finish = 250;
   else
      r.st.finish = finish;
   end
   if isempty (cues)
      r.st.cues = 23;
   else
      r.st.cues = cues;
   end
   if isempty (conditions)
      r.st.conds = -1;   
   else
      r.st.conds = conditions;
   end
   if isempty (blocks)
      r.st.blocks = -1;
   else
      r.st.blocks = blocks;
   end
   if isempty (responses)
      r.st.resps = -1;
   else
      r.st.resps = responses;
   end
   if isempty (trials)
      r.st.trials = -1;
   else
      r.st.trials = trials;
   end
   r.size = 1;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1 % the argument is a file name
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   f = openfile (label, 'r');
   
   % some useful variables
   total_rules     = 0;
   rules_per_group = 0;
   group           = 0;
   old_label       ='';
   num_of_fields   = 0;
   
   % read the first line and see what fields are included
   head = deblank (fgetl (f));
   while ~isempty (head)
      num_of_fields = num_of_fields + 1;
      [field{num_of_fields}, head] = strtok (head);
   end
   if ~strcmp (field {1}, 'label')
      error ('''RULES'' ERROR: first field in groups file must be ''label''');
   end
   
   % read the rest of the file, line by line
   while 1
      % read one line
      lin = fgetl (f);
      
      % check for feof
      if lin == -1
         break
      end
      
      oneline = deblank (lin);
      [value, oneline] = strtok (oneline);
      
      % if (label) then there is something to process; 
      % otherwise oneline is just an empty line,
      % in which case skip the whole thing and get the next
      if value & value(1)~='%'
         
         % if the current rule label is different from the previous
         % then this is a new group (also applies to the very first rule)
         
         total_rules = total_rules + 1;
         r.st(total_rules).label = value;
         if ~strcmp (value, old_label)
            group = group + 1;
            r.size (group) = 0;
         end
         r.size (group) = r.size (group) + 1;
         old_label = value;
         
         % every field in the current rule is at first set to a default value
         r.st(total_rules).start  = -50;
         r.st(total_rules).finish = 250;
         r.st(total_rules).cues   =  23;
         r.st(total_rules).conds  =  -1;
         r.st(total_rules).blocks =  -1;
         r.st(total_rules).resps  =   0;
         r.st(total_rules).trials =  -1;
         
         % loop through the fields actually present in the file
         % except for 'label', which is compulsory
         
         for i=2:num_of_fields
            [value, oneline] = strtok (oneline);
            if isempty (value)
               error (['''RULES'' ERROR: too few values in rule ', total_rules]);
            end
            v = eval (value); %, error ('''RULES'' ERROR: syntax error in rule'));
            % try to assign the value to a 'rules' field
            switch (field{i}(1:3))
            case 'tri'
               r.st(total_rules).trials = v;
            case 'con'
               r.st(total_rules).conds = v;
            case 'blo'
               r.st(total_rules).blocks = v;
            case 'res'
               r.st(total_rules).resps = v;
            case 'cue'
               r.st(total_rules).cues = v;
            case 'sta'
               if length (v) > 1
                  error ('''RULES'' ERROR: START field must be a scalar');
               end
               r.st(total_rules).start = v;
            case 'fin'
               if length (v) > 1
                  error ('''RULES'' ERROR: FINISH field must be a scalar');
               end
               r.st(total_rules).finish = v;
            otherwise
               error ('''RULES'' ERROR: bad field in header');
            end
         end % resume here after having processed all fields in oneline
      end % resume here if oneline is empty
      if value & value(1) ~= '%'
         [value, oneline] = strtok (oneline);
         if value & value(1) ~= '%'
            t = num2str (total_rules);
            error (['''RULES'' ERROR: too many values in rule ', t]);
         end
      end
   end
   fclose (f);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   error ('''RULES'' ERROR: wrong number of arguments');
end

%finish up
r = class (r, 'rules');
