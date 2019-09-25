% plot the clusters for a decoding event
%event = 38;
event = 11
i = [3 1 7]
aname = 'frank';
basedir = '/data/mkarlsso/';

clist = [ 8     6     1     1;
     8     6     6     4;
     8     6     5     1;
     8     6    10     1;
     8     6    11     1;
     8     6    29     1;
     8     6    21     5;
     8     6     8     3;
     8     6     8     4;
     8     6    10     3;
     8     6     1     4];

clist(:,2) = 2 ;
plotclusterlist('frank', basedir, clist, [4 3]);
orient tall

figure
clist(:,2) = 6 ;
plotclusterlist('frank', basedir, clist, [4 3]);
orient tall

event = 163
i = [2 1 2]
aname = 'bond';
basedir = '/data/mkarlsso/';
clist =     [ 4     2    11     4 ;
	      4     2    12     1 ;
	      4     2    14     5 ;
	      4     2    13     1 ;
	      4     2    14     3 ;
	      4     2    10     1 ;
	      4     2    12     3 ;
	      4     2     5     1 ;
	      4     2    18     1 ;
	      4     2     1     1 ;
	      4     2     5     2 ;
	      4     2    19     2 ;
	      4     2     1    10 ;
	      4     2    17     1 ;
	      4     2    19     1 ;
	      4     2    29     4 ;
	      4     2     2     4 ;
	      4     2    11     5 ]


plotclusterlist('bond', basedir, clist, [6 3]);
orient tall

figure
clist(:,2) = 3;
plotclusterlist('bond', basedir, clist, [6 3]);
orient tall     


E1 in E2 #2
aname = 'bond';
basedir = '/data/mkarlsso/';
clist =     [ 4     6    12     1;
     4     6    11     4;
     4     6    13     1;
     4     6    14     3;
     4     6    10     1;
     4     6    12     3;
     4     6    18     1;
     4     6     5     1;
     4     6    14     1;
     4     6     5     2;
     4     6    17     1;
     4     6    19     1;
     4     6    19     2;
     4     6    14     5];



clist(:,2) = 2;
plotclusterlist('bond', basedir, clist, [5 3]);
orient tall

figure
clist(:,2) = 6;
plotclusterlist('bond', basedir, clist, [5 3]);
orient tall     

