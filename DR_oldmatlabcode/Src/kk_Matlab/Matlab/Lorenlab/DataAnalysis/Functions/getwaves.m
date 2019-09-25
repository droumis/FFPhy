function [outwaves,outts] = getwaves(animal,indices)

N = size(indices,1);

for n = 1:N
   day = indices(n,1);
   adir = animaldirdata(animal,day);
   if isempty(adir)
      error('Animal/day not found.');
   end
   d = dir(fullfile(adir,'*-*'));
   d(~[d.isdir]) = [];

   dir_tets = strvcat(d.name);
   dir_tets( ~ismember(dir_tets(:,1), sprintf('%d',[0:9])), :) = [];
   dir_tets = str2num(dir_tets(:,1:2));

   [tmp,adir_tail] = fileparts(adir);
   tetdir = d( find(dir_tets == indices(n,3)) ).name;
   load(fullfile(adir,tetdir,strcat(adir_tail,'-',tetdir(1:2),'.mat')));

   animdef = animaldef(animal);
   spikes = loaddatastruct(animdef{2},animdef{3},'spikes',day);
   ss = spikes{indices(n,1)}{indices(n,2)}{indices(n,3)}{indices(n,4)};
   ts = round(ss.data(:,1)*10000);

   inds = find(ismember(timestamps,ts));

   outwaves{n} = waves(:,:,inds);
   outts{n} = timestamps(inds);
end

if N == 1
   outwaves = outwaves{1};
   outts =  outts{1};
end

