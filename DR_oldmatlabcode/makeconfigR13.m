%for 18-tetrode drive, loaded continuously, missing 7, 22, 30
%8, 15, and 23 are tetrodes loaded across chips
%last three electrodes not displayed on EEG, since isn't room, leaving 24
%load anamap
map1 = [1 1 1 1 2 2 2 2 3 3 3 3 0 0 4 4 4 4 5 5 5 5 6 6 6 6 7]';
map2 = [8 8 8 8 9 9 9 9 10 10 10 10 11 11 11 11 12 12 12 12 13 13 13 13 14 14
14]';
%map3 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 23]';
map4 = [23 23 23 23 24 24 24 24 25 25 25 25 26 26 26 26 27 27 27 27 28 28 28 29 29 29 29]';
anamap = [map1 map2 map4];
GenerateDSPmap(anamap, 'specialtetrodes',[15],'master','ike','slave1','mike','slave2','spike','rt','tyke', 'datadir','/data/ana/')
%GenerateDSPmap(anamap,'specialtetrodes',[],'eegexclude',[1 2 3 5 10 11 20 23 25 26 29],'master','drizzle','slave1','rain','slave2','mist','rt','fog','datadir','/data/mkarlsso/')
