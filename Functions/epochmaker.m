% returns full epochfilter for given epochstring handle
function epochfilt = epochmaker(varargin)

vars = {varargin{1:2:end}};
vals = {varargin{2:2:end}};
epochfilt = '';
for ivar == 1:length(vars)
    if ~strfind(ivar, 'num') %deal with things like numca1tets and numgammatets etc
        epochfilt = [epochfilt sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',ivar, ival);
        eval(sprintf('%s=%s;',vars{i},mat2str(vals{i}))); %assign in var to val in workspcae
    else
    end
    
end

end


%     case 'runlinear'  % any run epoch in linear environment w/ minimum detecting tets
%         epochfilter = ['((($num_gammaftets >= 2) && ($numca1tets >= 3)) && (isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB'')) || (isequal($environment, ''LinearA'')) || (isequal($environment, ''LinearB''))  )'];
%     case 'runlinear_inclusive'
%         epochfilter = ['((isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB'')) || (isequal($environment, ''LinearA'')) || (isequal($environment, ''LinearB''))  )'];       
%     case 'runlinear_rip'
%         epochfilter = ['($numca1tets >= 3) && ((isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB'')) || (isequal($environment, ''LinearA'')) || (isequal($environment, ''LinearB''))  )'];              
%     case 'runW_rip'
%         epochfilter = ['($numca1tets >= 3) && ((isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB''))  )'];              
%     case 'runW'       % any run epoch in a W-track environment w/ minimum detecting tets
%         epochfilter = ['((($num_gammaftets >= 2) && ($numca1tets >= 3)) && (isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB''))   )'];
%     case 'run'
%         epochfilter = ['(($numca1tets >= 3) && ($num_gammaftets >= 2) && (  isequal($type, ''run'') )   )'];
%     case 'sleep'
%         epochfilter = ['(($numca1tets >= 3) && ($num_gammaftets >= 2) && (  isequal($type, ''sleep'') )   )'];
%     case 'all'
%         epochfilter = ['(($numca1tets >= 3) && ($num_gammaftets >= 2) && (isequal($type, ''run'') || isequal($type, ''sleep'')))'];
%     case 'all2'
%         epochfilter = ['(isequal($type, ''run'') || isequal($type, ''sleep''))'];
%     case 'all_rip'
%         epochfilter = ['(($numca1tets >= 3) && (isequal($type, ''run'') || isequal($type, ''sleep'')))'];
%     case 'sleep_inclusive'
%         epochfilter = ['(($numca1tets >= 3) && (  isequal($type, ''sleep'') )   )'];        
%     case 'runlinear_all+sleep_inclusive'
%         epochfilter = [ ['((($numca1tets >= 3) && (  isequal($type, ''sleep'') )   )) || '] ...
%                         ['(((isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB'')) || (isequal($environment, ''LinearA'')) || (isequal($environment, ''LinearB''))  ))']];           
%     case 'runW_rip+sleep_inclusive'
%         epochfilter = [ ['((($numca1tets >= 3) && (  isequal($type, ''sleep'') )   )) || '] ...
%                         ['(((isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB''))  ))']];                  
%     case 'runlinear_all'  % any run epoch in linear environment
%         epochfilter = ['((isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB'')) || (isequal($environment, ''LinearA'')) || (isequal($environment, ''LinearB''))  )'];
%     case 'runW_rip'  % W-track epochs w/ sufficient ca1 detecting tets
%         epochfilter = ['( ($numca1tets >= 3) && ((isequal($environment, ''WTrackA'')) || (isequal($environment, ''WTrackB'')) || (isequal($environment, ''TrackA'')) || (isequal($environment, ''TrackB'')) ) )'];
%     case 'openfield'  % Open Field
%         epochfilter = ['( ($numca1tets >= 3) && ((isequal($environment, ''OpenFieldA'')) || (isequal($environment, ''OpenFieldB''))  ) )'];
%     case 'everyep'
%         epochfilter = [' ( isequal($type, ''run'') || isequal($type, ''sleep'')  )'] ;
%     case 'run_all'
%         epochfilter = ['(isequal($type, ''run''))'];
%     case 'sleep_all'
%         epochfilter = ['(isequal($type, ''sleep''))'];
%     case