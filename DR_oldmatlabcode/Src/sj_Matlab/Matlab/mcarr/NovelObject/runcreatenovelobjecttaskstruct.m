%Create Novel Object Task Struct
createnovelobjecttaskstruct('/data6/monster/Cor/','Cor',[10 2; 10 4; 10 6; 10 8; 11 2; 11 4; 11 6; 11 8], [82 51])
createnovelobjecttaskstruct('/data6/monster/Cyc/','Cyc',[1 2; 1 4; 1 6; 2 2; 2 4; 2 6; 3 2; 3 4; 3 6; 4 2; 4 4; 4 6; 5 2; 5 4; 5 6; 6 2; 6 4; 6 6], [66 52])
createnovelobjecttaskstruct('/data6/monster/Faf/','Faf',[1 2; 1 4; 1 6; 2 2; 2 4; 2 6; 3 2; 3 4; 3 6; 4 2; 4 4; 4 6; 5 2; 5 4; 5 6; 6 2; 6 4; 6 6], [68 57])
createnovelobjecttaskstruct('/data6/monster/Gre/','Gre',[1 2; 1 4; 1 6; 2 2; 2 4; 2 6; 3 2; 3 4; 3 6; 4 2; 4 4; 4 6; 5 2; 5 4; 5 6; 6 2; 6 4; 6 6; 7 2; 7 4; 7 6], [61 52])
createnovelobjecttaskstruct('/data6/monster/Dun/','Dun',[2 2; 2 4; 2 6; 2 10; 2 12; 2 14;4 2; 4 4; 4 6; 5 2; 5 4; 5 6; 6 2; 6 4; 6 6; 7 2; 7 4; 7 6], [57 56])
createnovelobjecttaskstruct('/data6/monster/God/','God',[1 2; 1 4; 1 6], [65 48.5])
createnovelobjecttaskstruct('/data6/monster/Cml/','Cml',[1 2; 1 4; 1 6; 1 8],[113.5 73.5])
createnovelobjecttaskstruct('/data6/monster/Cml/','Cml',[2 2; 2 4; 2 6; 2 8],[117.5 74])
createnovelobjecttaskstruct('/data6/monster/Cml/','Cml',[3 2; 3 4; 3 6; 3 8],[119.5 80])
createnovelobjecttaskstruct('/data6/monster/Cml/','Cml',[4 2; 4 4; 4 6; 4 8],[118.5 81.5])
createnovelobjecttaskstruct('/data6/monster/Nic/','Nic',[1 2; 1 4; 1 6; 1 8],[109 76.5])
createnovelobjecttaskstruct('/data6/monster/Nic/','Nic',[2 2; 2 4; 2 6; 2 8],[110 77.5])
createnovelobjecttaskstruct('/data6/monster/Nic/','Nic',[3 2; 3 4; 3 6; 3 8],[106 81.5])   
createnovelobjecttaskstruct('/data6/monster/Nic/','Nic',[4 2; 4 4; 4 6; 4 8],[108 83])
createnovelobjecttaskstruct('/data6/monster/Nic/','Nic',[5 2; 5 4; 5 6; 5 8],[106.5 93])
createnovelobjecttaskstruct('/data6/monster/Nic/','Nic',[6 2; 6 4; 6 6; 6 8],[115.5 91.5])
createnovelobjecttaskstruct('/data6/monster/Nic/','Nic',[7 2; 7 4; 7 6; 7 8],[115 88.5])
createnovelobjecttaskstruct('/data6/monster/Nic/','Nic',[8 2; 8 4; 8 6; 8 8],[115.5 79.5])

%Add Description of session: familiar, novel, or super novel
%Cor
addtaskinfo('/data6/monster/Cor/','Cor',10:11,[2 6],'session','familiar')
addtaskinfo('/data6/monster/Cor/','Cor',10:11,[4 8],'session','novel')
addtaskinfo('/data6/monster/Cor/','Cor',10:11,[1 3 5 7 9],'type','sleep')
addtaskinfo('/data6/monster/Cor/','Cor',10:11,1,'runbefore','baseline');
addtaskinfo('/data6/monster/Cor/','Cor',10:11,[3 7],'runbefore','familiar');
addtaskinfo('/data6/monster/Cor/','Cor',10:11,[5 9],'runbefore','novel');

%Cyc
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,2,'session','familiar')
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,4,'session','novel')
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,6,'session','supernovel')
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,[1 3 5 7],'type','sleep')
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,1,'runbefore','baseline');
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,3,'runbefore','familiar');
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,5,'runbefore','novel');
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,7,'runbefore','supernovel');

%Dun
addtaskinfo('/data6/monster/Dun/','Dun',[2 4 5 6 7],2,'session','familiar')
addtaskinfo('/data6/monster/Dun/','Dun',[2 4 5 6 7],4,'session','novel')
addtaskinfo('/data6/monster/Dun/','Dun',[2 4 5 6 7],6,'session','supernovel')
addtaskinfo('/data6/monster/Dun/','Dun',2,10,'session','familiar')
addtaskinfo('/data6/monster/Dun/','Dun',2,12,'session','novel')
addtaskinfo('/data6/monster/Dun/','Dun',2,14,'session','supernovel')
addtaskinfo('/data6/monster/Dun/','Dun',[2 4 5 6 7],[1 3 5 7],'type','sleep')
addtaskinfo('/data6/monster/Dun/','Dun',[2 4 5 6 7],1,'runbefore','baseline');
addtaskinfo('/data6/monster/Dun/','Dun',[2 4 5 6 7],3,'runbefore','familiar');
addtaskinfo('/data6/monster/Dun/','Dun',[2 4 5 6 7],5,'runbefore','novel');
addtaskinfo('/data6/monster/Dun/','Dun',[2 4 5 6 7],7,'runbefore','supernovel');
addtaskinfo('/data6/monster/Dun/','Dun',2,[9 11 13 15],'type','sleep')
addtaskinfo('/data6/monster/Dun/','Dun',2,9,'runbefore','baseline');
addtaskinfo('/data6/monster/Dun/','Dun',2,11,'runbefore','familiar');
addtaskinfo('/data6/monster/Dun/','Dun',2,13,'runbefore','novel');
addtaskinfo('/data6/monster/Dun/','Dun',2,15,'runbefore','supernovel');

%Faf
addtaskinfo('/data6/monster/Faf/','Faf',1:6,2,'session','familiar')
addtaskinfo('/data6/monster/Faf/','Faf',1:6,4,'session','novel')
addtaskinfo('/data6/monster/Faf/','Faf',1:6,6,'session','supernovel')
addtaskinfo('/data6/monster/Faf/','Faf',1:6,[1 3 5 7],'type','sleep')
addtaskinfo('/data6/monster/Faf/','Faf',1:6,1,'runbefore','baseline');
addtaskinfo('/data6/monster/Faf/','Faf',1:6,3,'runbefore','familiar');
addtaskinfo('/data6/monster/Faf/','Faf',1:6,5,'runbefore','novel');
addtaskinfo('/data6/monster/Faf/','Faf',1:6,7,'runbefore','supernovel');

%God
addtaskinfo('/data6/monster/God/','God',1,2,'session','familiar')
addtaskinfo('/data6/monster/God/','God',1,4,'session','novel')
addtaskinfo('/data6/monster/God/','God',1,6,'session','supernovel')
addtaskinfo('/data6/monster/God/','God',1,[1 3 5 7],'type','sleep')
addtaskinfo('/data6/monster/God/','God',1,1,'runbefore','baseline');
addtaskinfo('/data6/monster/God/','God',1,3,'runbefore','familiar');
addtaskinfo('/data6/monster/God/','God',1,5,'runbefore','novel');
addtaskinfo('/data6/monster/God/','God',1,7,'runbefore','supernovel');

%Gre
addtaskinfo('/data6/monster/Gre/','Gre',1:7,2,'session','familiar')
addtaskinfo('/data6/monster/Gre/','Gre',1:7,4,'session','novel')
addtaskinfo('/data6/monster/Gre/','Gre',1:7,6,'session','supernovel')
addtaskinfo('/data6/monster/Gre/','Gre',1:7,[1 3 5 7],'type','sleep')
addtaskinfo('/data6/monster/Gre/','Gre',1:7,1,'runbefore','baseline');
addtaskinfo('/data6/monster/Gre/','Gre',1:7,3,'runbefore','familiar');
addtaskinfo('/data6/monster/Gre/','Gre',1:7,5,'runbefore','novel');
addtaskinfo('/data6/monster/Gre/','Gre',1:7,7,'runbefore','supernovel');

%Cml
addtaskinfo('/data6/monster/Cml/','Cml',1:4,2,'session','familiar_repeat')
addtaskinfo('/data6/monster/Cml/','Cml',1:4,4,'session','familiar')
addtaskinfo('/data6/monster/Cml/','Cml',1:4,6,'session','novel')
addtaskinfo('/data6/monster/Cml/','Cml',1:4,8,'session','supernovel')
addtaskinfo('/data6/monster/Cml/','Cml',1:4,[1 3 5 7 9],'type','sleep')
addtaskinfo('/data6/monster/Cml/','Cml',1:4,[2 4 6 8],'type','run')
addtaskinfo('/data6/monster/Cml/','Cml',1:4,1,'runbefore','baseline');
addtaskinfo('/data6/monster/Cml/','Cml',1:4,5,'runbefore','familiar');
addtaskinfo('/data6/monster/Cml/','Cml',1:4,7,'runbefore','novel');
addtaskinfo('/data6/monster/Cml/','Cml',1:4,9,'runbefore','supernovel');

%Nico
addtaskinfo('/data6/monster/Nic/','Nic',1:8,2,'session','familiar_repeat')
addtaskinfo('/data6/monster/Nic/','Nic',1:8,4,'session','familiar')
addtaskinfo('/data6/monster/Nic/','Nic',1:8,6,'session','novel')
addtaskinfo('/data6/monster/Nic/','Nic',1:8,8,'session','supernovel')
addtaskinfo('/data6/monster/Nic/','Nic',1:8,[1 3 5 7 9],'type','sleep')
addtaskinfo('/data6/monster/Nic/','Nic',1:8,[2 4 6 8],'type','run')
addtaskinfo('/data6/monster/Nic/','Nic',1:8,1,'runbefore','baseline');
addtaskinfo('/data6/monster/Nic/','Nic',1:8,5,'runbefore','familiar');
addtaskinfo('/data6/monster/Nic/','Nic',1:8,7,'runbefore','novel');
addtaskinfo('/data6/monster/Nic/','Nic',1:8,9,'runbefore','supernovel');

%Add object information for each animal
%Cor
addtaskinfo('/data6/monster/Cor/','Cor',10,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cor/','Cor',10,4,'objects',[0 0 1 1],'novel',[0 1 1 0],'familiarsession',2)
addtaskinfo('/data6/monster/Cor/','Cor',10,6,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cor/','Cor',10,8,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',6)
addtaskinfo('/data6/monster/Cor/','Cor',11,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cor/','Cor',11,4,'objects',[0 1 1 0],'novel',[0 0 1 1],'familiarsession',2)
addtaskinfo('/data6/monster/Cor/','Cor',11,6,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cor/','Cor',11,8,'objects',[1 0 0 1],'novel',[1 1 0 0],'familiarsession',6)

%Cyc
addtaskinfo('/data6/monster/Cyc/','Cyc',1:6,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cyc/','Cyc',1,4,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',2)
addtaskinfo('/data6/monster/Cyc/','Cyc',1,6,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Cyc/','Cyc',2,4,'objects',[0 0 1 1],'novel',[0 1 1 0],'familiarsession',2)
addtaskinfo('/data6/monster/Cyc/','Cyc',2,6,'objects',[1 0 0 1],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Cyc/','Cyc',3,4,'objects',[0 1 1 0],'novel',[0 0 1 1],'familiarsession',2)
addtaskinfo('/data6/monster/Cyc/','Cyc',3,6,'objects',[1 1 0 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Cyc/','Cyc',4,4,'objects',[1 0 0 1],'novel',[1 1 0 0],'familiarsession',2)
addtaskinfo('/data6/monster/Cyc/','Cyc',4,6,'objects',[1 0 0 1],'novel',[1 0 0 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Cyc/','Cyc',5,4,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',2)
addtaskinfo('/data6/monster/Cyc/','Cyc',5,6,'objects',[1 0 0 1],'novel',[0 1 0 1],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Cyc/','Cyc',6,4,'objects',[1 0 0 1],'novel',[1 1 0 0],'familiarsession',2)
addtaskinfo('/data6/monster/Cyc/','Cyc',6,6,'objects',[0 0 1 1],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)

%Faf
addtaskinfo('/data6/monster/Faf/','Faf',1:6,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Faf/','Faf',1,4,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',2)
addtaskinfo('/data6/monster/Faf/','Faf',2,4,'objects',[0 0 1 1],'novel',[0 1 1 0],'familiarsession',2)
addtaskinfo('/data6/monster/Faf/','Faf',3,4,'objects',[0 1 1 0],'novel',[0 0 1 1],'familiarsession',2)
addtaskinfo('/data6/monster/Faf/','Faf',4,4,'objects',[1 0 0 1],'novel',[1 1 0 0],'familiarsession',2)
addtaskinfo('/data6/monster/Faf/','Faf',5,4,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',2)
addtaskinfo('/data6/monster/Faf/','Faf',6,4,'objects',[0 0 1 1],'novel',[0 1 1 0],'familiarsession',2)
addtaskinfo('/data6/monster/Faf/','Faf',1,6,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Faf/','Faf',2,6,'objects',[1 0 0 1],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Faf/','Faf',3,6,'objects',[1 1 0 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Faf/','Faf',4,6,'objects',[0 0 1 1],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Faf/','Faf',5,6,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Faf/','Faf',6,6,'objects',[1 0 0 1],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)

%Gre
addtaskinfo('/data6/monster/Gre/','Gre',1:7,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Gre/','Gre',1,4,'objects',[0 1 1 0],'novel',[0 0 1 1],'familiarsession',2)
addtaskinfo('/data6/monster/Gre/','Gre',1,6,'objects',[1 1 0 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Gre/','Gre',2,4,'objects',[1 0 0 1],'novel',[1 1 0 0],'familiarsession',2)
addtaskinfo('/data6/monster/Gre/','Gre',2,6,'objects',[0 0 1 1],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Gre/','Gre',3,4,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',2)
addtaskinfo('/data6/monster/Gre/','Gre',3,6,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Gre/','Gre',4,4,'objects',[0 0 1 1],'novel',[0 1 1 0],'familiarsession',2)
addtaskinfo('/data6/monster/Gre/','Gre',4,6,'objects',[1 0 0 1],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Gre/','Gre',5,4,'objects',[0 1 1 0],'novel',[0 0 1 1],'familiarsession',2)
addtaskinfo('/data6/monster/Gre/','Gre',5,6,'objects',[1 1 0 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Gre/','Gre',6,4,'objects',[1 0 0 1],'novel',[1 1 0 0],'familiarsession',2)
addtaskinfo('/data6/monster/Gre/','Gre',6,6,'objects',[0 0 1 1],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Gre/','Gre',7,4,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',2)
addtaskinfo('/data6/monster/Gre/','Gre',7,6,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)

%Dun
addtaskinfo('/data6/monster/Dun/','Dun',2,2,'objects',[1 0 1 0],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Dun/','Dun',2,4,'objects',[0 1 1 0],'novel',[1 1 0 0],'familiarsession',2)
addtaskinfo('/data6/monster/Dun/','Dun',2,6,'objects',[0 0 1 1],'novel',[0 1 0 1],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Dun/','Dun',2,10,'objects',[1 0 1 0],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Dun/','Dun',2,12,'objects',[1 1 0 0],'novel',[0 1 1 0],'familiarsession',10)
addtaskinfo('/data6/monster/Dun/','Dun',2,14,'objects',[1 0 0 1],'novel',[0 1 0 1],'familiarsession',10,'novelsession',12)
addtaskinfo('/data6/monster/Dun/','Dun',4,2,'objects',[1 0 1 0],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Dun/','Dun',4,4,'objects',[1 1 0 0],'novel',[0 1 1 0],'familiarsession',2)
addtaskinfo('/data6/monster/Dun/','Dun',4,6,'objects',[1 0 0 1],'novel',[0 1 0 1],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Dun/','Dun',5:7,2,'objects',[1 0 1 0],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Dun/','Dun',5,4,'objects',[0 1 1 0],'novel',[1 1 0 0],'familiarsession',2)
addtaskinfo('/data6/monster/Dun/','Dun',6,4,'objects',[1 0 0 1],'novel',[0 0 1 1],'familiarsession',2)
addtaskinfo('/data6/monster/Dun/','Dun',7,4,'objects',[0 0 1 1],'novel',[1 0 0 1],'familiarsession',2)
addtaskinfo('/data6/monster/Dun/','Dun',5,6,'objects',[0 0 1 1],'novel',[0 1 0 1],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Dun/','Dun',6,6,'objects',[1 1 0 0],'novel',[0 1 0 1],'familiarsession',2,'novelsession',4)
addtaskinfo('/data6/monster/Dun/','Dun',7,6,'objects',[0 1 1 0],'novel',[0 1 0 1],'familiarsession',2,'novelsession',4)

%God
addtaskinfo('/data6/monster/God/','God',1,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/God/','God',1,4,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',2)
addtaskinfo('/data6/monster/God/','God',1,6,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',2,'novelsession',4)

%Cml
addtaskinfo('/data6/monster/Cml/','Cml',1,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cml/','Cml',1,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cml/','Cml',1,6,'objects',[0 0 1 1],'novel',[0 1 1 0],'familiarsession',4)
addtaskinfo('/data6/monster/Cml/','Cml',1,8,'objects',[1 0 0 1],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Cml/','Cml',2,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cml/','Cml',2,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cml/','Cml',2,6,'objects',[1 0 0 1],'novel',[1 1 0 0],'familiarsession',4)
addtaskinfo('/data6/monster/Cml/','Cml',2,8,'objects',[0 0 1 1],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Cml/','Cml',3,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cml/','Cml',3,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cml/','Cml',3,6,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',4)
addtaskinfo('/data6/monster/Cml/','Cml',3,8,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Cml/','Cml',4,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cml/','Cml',4,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Cml/','Cml',4,6,'objects',[0 1 1 0],'novel',[0 0 1 1],'familiarsession',4)
addtaskinfo('/data6/monster/Cml/','Cml',4,8,'objects',[1 1 0 0],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)

%Nico
addtaskinfo('/data6/monster/Nic/','Nic',1,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',1,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',1,6,'objects',[0 0 1 1],'novel',[0 1 1 0],'familiarsession',4)
addtaskinfo('/data6/monster/Nic/','Nic',1,8,'objects',[1 0 0 1],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Nic/','Nic',2,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',2,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',2,6,'objects',[1 0 0 1],'novel',[1 0 0 1],'familiarsession',4)
addtaskinfo('/data6/monster/Nic/','Nic',2,8,'objects',[0 0 1 1],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Nic/','Nic',3,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',3,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',3,6,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',4)
addtaskinfo('/data6/monster/Nic/','Nic',3,8,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Nic/','Nic',4,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',4,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',4,6,'objects',[0 1 1 0],'novel',[0 0 1 1],'familiarsession',4)
addtaskinfo('/data6/monster/Nic/','Nic',4,8,'objects',[1 1 0 0],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Nic/','Nic',5,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',5,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',5,6,'objects',[0 0 1 1],'novel',[0 1 1 0],'familiarsession',4)
addtaskinfo('/data6/monster/Nic/','Nic',5,8,'objects',[1 0 0 1],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Nic/','Nic',6,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',6,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',6,6,'objects',[1 0 0 1],'novel',[1 0 0 1],'familiarsession',4)
addtaskinfo('/data6/monster/Nic/','Nic',6,8,'objects',[0 0 1 1],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Nic/','Nic',7,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',7,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',7,6,'objects',[1 1 0 0],'novel',[1 0 0 1],'familiarsession',4)
addtaskinfo('/data6/monster/Nic/','Nic',7,8,'objects',[0 1 1 0],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)
addtaskinfo('/data6/monster/Nic/','Nic',8,2,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',8,4,'objects',[0 1 0 1],'novel',[0 0 0 0])
addtaskinfo('/data6/monster/Nic/','Nic',8,6,'objects',[0 1 1 0],'novel',[0 0 1 1],'familiarsession',4)
addtaskinfo('/data6/monster/Nic/','Nic',8,8,'objects',[1 1 0 0],'novel',[1 0 1 0],'familiarsession',4,'novelsession',6)

%% NOW MAKE CELLINFOSTRUCT AND TETINFOSTRUCT

createcellinfostruct('/data6/monster/Cor/','Cor')
createcellinfostruct('/data6/monster/Cyc/','Cyc')
createcellinfostruct('/data6/monster/Faf/','Faf')
createcellinfostruct('/data6/monster/Gre/','Gre')
createcellinfostruct('/data6/monster/Dun/','Dun')
createcellinfostruct('/data6/monster/God/','God')
createcellinfostruct('/data6/monster/Cml/','Cml')
createcellinfostruct('/data6/monster/Nic/','Nic')

createtetinfostruct('/data6/monster/Cor/','Cor')
createtetinfostruct('/data6/monster/Cyc/','Cyc')
createtetinfostruct('/data6/monster/Faf/','Faf')
createtetinfostruct('/data6/monster/Gre/','Gre')
createtetinfostruct('/data6/monster/Dun/','Dun')
createtetinfostruct('/data6/monster/God/','God')
createtetinfostruct('/data6/monster/Cml/','Cml')
createtetinfostruct('/data6/monster/Nic/','Nic')

addtetrodelocation('/data6/monster/Cor/','Cor',[16 17 18 19 20 21 22 23],'CA1')
addtetrodelocation('/data6/monster/Cor/','Cor',[14 15 25],'CA3')
addtetrodelocation('/data6/monster/Cor/','Cor',24,'Reference')

addtetrodelocation('/data6/monster/Cyc/','Cyc',[1 2 3 8 9 10 11 12],'CA1')
addtetrodelocation('/data6/monster/Cyc/','Cyc',[4 5],'CA3')
addtetrodelocation('/data6/monster/Cyc/','Cyc',7,'Reference')

addtetrodelocation('/data6/monster/Dun/','Dun',[2 4 6 12],'CA1')
addtetrodelocation('/data6/monster/Dun/','Dun',[5 7 9 10 11 13],'CA3')
addtetrodelocation('/data6/monster/Dun/','Dun',3,'Reference')

addtetrodelocation('/data6/monster/Faf/','Faf',[1 2 4 6 12 13],'CA1')
addtetrodelocation('/data6/monster/Faf/','Faf',[8 9 10 14],'CA3')
addtetrodelocation('/data6/monster/Faf/','Faf',11,'Reference')

addtetrodelocation('/data6/monster/God/','God',[2 9 10 11 12 13 14 15],'CA1')
addtetrodelocation('/data6/monster/God/','God',[1 4 6 7 8],'CA3')
addtetrodelocation('/data6/monster/God/','God',3,'Reference')

addtetrodelocation('/data6/monster/Gre/','Gre',[1 2 3 4 5 6 7 9 13],'CA1')
addtetrodelocation('/data6/monster/Gre/','Gre',12,'CA2')
addtetrodelocation('/data6/monster/Gre/','Gre',[10 11 14],'CA3')
addtetrodelocation('/data6/monster/Gre/','Gre',8,'Reference')

addtetrodelocation('/data6/monster/Cml/','Cml',[3 6 9 12 13 14],'CA1')
addtetrodelocation('/data6/monster/Cml/','Cml',[2 4 5 ],'CA3')
addtetrodelocation('/data6/monster/Cml/','Cml',1,'Reference')

addtetrodelocation('/data6/monster/Nic/','Nic',[1 3 6 7 8 9 11 13],'CA1')
addtetrodelocation('/data6/monster/Nic/','Nic',[4 12],'CA3')
addtetrodelocation('/data6/monster/Nic/','Nic',5,'Reference')

