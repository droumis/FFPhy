function trajdayprocess(directoryname,fileprefix,days)
%TRAJDAYPROCESS(directoryname,fileprefix,days, options)
%
% Runs getbehavestate for all run epochs in each day and saves the data in
%'linpos.traj' in the directoryname folder.  See GETBEHAVESTATE for more
%details
%
%directoryname - example 'data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile
%
%days -          a vector of experiment day numbers 
%

lowercasethree = '';

days = days(:)';

for day = days
   
   dsz = '';
   if (day < 10)
      dsz = '0';
   end

   eval(['load ',directoryname,fileprefix,'linpos',dsz,num2str(day), '.mat']);
   eval(['load ',directoryname,fileprefix,'task',dsz,num2str(day), '.mat']);
   eval(['task = ',lowercasethree,'task;'])
   for i = 1:length(task{day})   
      if ((~isempty(task{day}{i})) && (strcmp(task{day}{i}.type,'run')) )
            
            disp(['Day ',num2str(day), ', Epoch ',num2str(i)])
            index = [day i];
            linpos{day}{i}.statematrix.traj = cell(1,6);
            for j=1:6
                [state, lindist] = getbehavestate(linpos,day,i,j);
                linpos{day}{i}.statematrix.traj{j} = state;
            end
      end
   end
   
   eval([lowercasethree,'linpos = linpos;']);
   eval(['save ',directoryname,fileprefix,'linpos',dsz,num2str(day),' ',lowercasethree,'linpos']);
end