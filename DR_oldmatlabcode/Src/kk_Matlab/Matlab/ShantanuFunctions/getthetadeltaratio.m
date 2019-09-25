function out = getthetadeltaratio(animaldir,animalprefix,epochs,tetrode)
% out = getthetadeltaratio(animaldir,animalprefix,epochs,tetrode)
% Produces a cell structure with the theta to delta ratio for each epoch
% tetrode
% EPOCHS - N by 2 matrix, columns are [day epoch]
% TETRODE - N by 1 matrix with tetrode which should be used for that day
% epoch pair

%Organize list of tetrodes to be used for each epoch
tet = [];
for j=1:size(tetrode,2)
    tet = [tet tetrode{j}];
end
tetrode = tet;


%Load data and compute theta to delta ratio
for i=1:size(epochs,1)
    theta = loadeegstruct(animaldir, animalprefix, 'theta', epochs(i,1), epochs(i,2), tetrode{i});
    delta = loadeegstruct(animaldir, animalprefix, 'delta', epochs(i,1), epochs(i,2), tetrode{i});
    time = geteegtimes(theta{epochs(i,1)}{epochs(i,2)}{tetrode{i}});
    
    theta = double(theta{epochs(i,1)}{epochs(i,2)}{tetrode{i}}.data(:,3));
    delta = double(delta{epochs(i,1)}{epochs(i,2)}{tetrode{i}}.data(:,3));
    
    %Convolve with a gaussian of 1 std and a width of 1 second
    g = gaussian(1,150);
    t = conv(theta,g);
    d = conv(delta,g);
    
    t = t(1:end-length(g)+1);
    d = d(1:end-length(g)+1);
    
    %Fill in any zeros in d with a very small number and compute theta to 
    %delta ratio
    d(d==0) = 0.00001;
    ratio = t./d;
    
    out{epochs(i,1)}{epochs(i,2)}.ratio = ratio';
    out{epochs(i,1)}{epochs(i,2)}.time = time;
end
