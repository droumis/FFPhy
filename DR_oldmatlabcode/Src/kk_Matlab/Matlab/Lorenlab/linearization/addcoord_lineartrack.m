function addcoord_lineartrack(animal)

[aname,adir,aprefx] = animaldef(animal);

pos = loaddatastruct(adir,aprefx,'pos');
task2 = loaddatastruct(adir,aprefx,'task');
epochs = evaluatefilter(task2,'strcmp($type,''run'')');

for e = 1:size(epochs,1)
  coord = getcoord_ntrack([],pos{epochs(e,1)}{epochs(e,2)}.data,[]);
  task2{epochs(e,1)}{epochs(e,2)}.linearcoord = coord;
  task{epochs(e,1)} = task2{epochs(e,1)};
  % save task
  filename = fullfile(adir,sprintf('%stask%02d.mat',aprefx,epochs(e,1)));
  save(filename,'task');
  clear task;
end

