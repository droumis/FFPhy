%% CMPERPIX
%% TETINFO
%% CELLINFO
%% UNIT CLASSIFICATION
%% TASKSTRUCT


    %% CHAPATI

% cmperpix
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati01','/datatmp/kkay/Cha','cha',1,'cmperpix',[.5 .5 .685 .685 .5 .5 .685 .685 .5])
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati02','/datatmp/kkay/Cha','cha',2,'cmperpix',[.5 .5 .685 .685 .5 .5 .685 .685 .5])
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati03','/datatmp/kkay/Cha','cha',3,'cmperpix',[.5 .685 .5 .685 .5])
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati04','/datatmp/kkay/Cha','cha',4,'cmperpix',[.5 .685 .5 .685 .5])
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati05','/datatmp/kkay/Cha','cha',5,'cmperpix',[.5 .685 .5 .685 .5])
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati06','/datatmp/kkay/Cha','cha',6,'cmperpix',[.5 .685 .5 .685 .5 .685 .5 .5])
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati07','/datatmp/kkay/Cha','cha',7,'cmperpix',[.5 .685 .5 .5 1.173])
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati08','/datatmp/kkay/Cha','cha',8,'cmperpix',[.5 .685 .5 .685 .5 1.173 .5])
cd('/datatmp/kkay/ChapatiData/')
kk_dayprocess_2('chapati09','/datatmp/kkay/Cha','cha',9,'cmperpix',[.5 .685 .5 .685 .5 1.173 .5])
cd('/datatmp/kkay/ChapatiData/')

% unit classification
unitclassifier_allprincipal('/datatmp/kkay/Cha/','cha',1:9,1:21)
typestruct_Chapati(1).type = 'inter';               % visual inspection 4.22.13
typestruct_Chapati(1).daytetcell = [ 1 12 1 ; 2 13 1 ; 3 6 1 ; ... % CA1     % visual inspection kk 4.22.13
                             4 9 8 ; 6 9 2 ;           ... % CA2
                             2 4 1 ];                   ... % CA3
unitclassifier_typestruct('/datatmp/kkay/Cha/','cha',1:9,1:21,typestruct_Chapati);

% taskstruct
% (already ran createtaskstruct and kk_updatetaskstruct_environment)
kk_updatetaskstruct_sleep('/datatmp/kkay/Cha/','cha',1:9)
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[1 3 1; 1 4 1; 1 7 2 ; 1 8 2])
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[2 3 3; 2 4 3; 2 7 4 ; 2 8 4])
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[3 2 5; 3 4 6])
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[4 2 7; 4 4 8])
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[5 2 9; 5 4 10])
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[6 2 11; 6 4 12; 6 6 13])
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[7 2 14])
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[8 2 16; 8 4 17; 8 6 2])
kk_updatetaskstruct_exposure('/datatmp/kkay/Cha/','cha',[9 2 18; 9 4 19; 9 6 3])


addgndtoeeg('cha',1,1:9,1:14)
addgndtoeeg('cha',2,1:9,1:14)
addgndtoeeg('cha',3,1:5,1:14)  % only went to 4 epochs..
addgndtoeeg('cha',4,1:5,1:14)
addgndtoeeg('cha',5,1:5,1:14)
addgndtoeeg('cha',6,1:8,1:14)
addgndtoeeg('cha',7,1:6,1:14)
addgndtoeeg('cha',8,1:7,1:14)
addgndtoeeg('cha',9,1:7,1:14)



    %% EGYPT

%cmperpix
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt01','/data12/mari/Egy','egy',1,'cmperpix',[.21 .5949 .21 .5949 .21 .21 .5949 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt02','/data12/mari/Egy','egy',2,'cmperpix',[.21 .5949 .21 .5949 .21 .5949 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt03','/data12/mari/Egy','egy',3,'cmperpix',[.21 .59 .21 .59 .21 .59 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt04','/data12/mari/Egy','egy',4,'cmperpix',[.21 .59 .21 .59 .21 .59 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt05','/data12/mari/Egy','egy',5,'cmperpix',[.21 .59 .21 .59 .21 .59 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt06','/data12/mari/Egy','egy',6,'cmperpix',[.21 .59 .21 .4812 .21 .21 .4812 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt07','/data12/mari/Egy','egy',7,'cmperpix',[.21 .4812 .21 .4812 .21 .4812 .21 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt08','/data12/mari/Egy','egy',8,'cmperpix',[.21 .4812 .21 .4812 .21 .4812 .21 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt09','/data12/mari/Egy','egy',9,'cmperpix',[.21 .479 .21 .479 .21 .479 .21 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt10','/data12/mari/Egy','egy',10,'cmperpix',[.21 .479 .21 .479 .21 .479 .21 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt11','/data12/mari/Egy','egy',11,'cmperpix',[.21 .479 .21 .479 .21 .479 .21])
cd('/data12/mari/EgyptData/')
kk_dayprocess_2('egypt12','/data12/mari/Egy','egy',12,'cmperpix',[.21 .479 .21 .479 .21 .479 .21 .21])

% unit classification
unitclassifier_allprincipal('/data12/mari/Egy/','egy',1:12,1:21)
typestruct_Egypt(1).type = 'inter';               % visual inspection kk 4.22.13
typestruct_Egypt(1).daytetcell = [ 2 15 4 ; 3 16 2 ; 4 16 4 ; 5 16 4 ; 6 16 2 ; 8 16 3 ; 9 16 7; ...  % CA1    visual inspection 4.22.13 kk
                             11 3 3 ; ...                                                       % CA2
                             1 10 1 ; 1 10 2 ; 2 10 1 ; 2 10 2 ; 3 14 4 ; 8 10 1 ];             % CA3  
typestruct_Egypt(2).type = 'triphasic';
typestruct_Egypt(2).daytetcell = [3 14 4 ; ...                                                                      % CA1   % visual inspection kk 4.22.13 kk
                            4 17 1 ; 5 17 3 ; 6 17 2 ; 7 17 3 ; 8 17 4 ; 9 17 4 ; 10 17 6 ; 11 17 3 ];    % CA2                                       
typestruct_Egypt(3).type = 'unknown';
typestruct_Egypt(3).daytetcell = [5 1 2];
unitclassifier_typestruct('/data12/mari/Egy/','egy',1:12,1:21,typestruct_Egypt);


createcellinfostruct('/data12/mari/Egy/','egy')
createtetinfostruct('/data12/mari/Egy/','egy')
addtetrodelocation_kk('/data12/mari/Egy/','egy',1,'CA3','CA3a')
addtetrodelocation_kk('/data12/mari/Egy/','egy',2,'CA1','CA1c')
addtetrodelocation_kk('/data12/mari/Egy/','egy',3,'CA2','CA2')
addtetrodelocation_kk('/data12/mari/Egy/','egy',6,'ctx','ctx-deep-overhpc')
addtetrodelocation_kk('/data12/mari/Egy/','egy',7,'CA3','CA3c')
addtetrodelocation_kk('/data12/mari/Egy/','egy',10,'CA3','CA3c')
addtetrodelocation_kk('/data12/mari/Egy/','egy',11,'cc','cc')
addtetrodelocation_kk('/data12/mari/Egy/','egy',13,'CA1','CA1c')
addtetrodelocation_kk('/data12/mari/Egy/','egy',14,'CA3','CA3b')
addtetrodelocation_kk('/data12/mari/Egy/','egy',15,'CA1','CA1c')
addtetrodelocation_kk('/data12/mari/Egy/','egy',16,'CA1','CA1-CA2')
addtetrodelocation_kk('/data12/mari/Egy/','egy',17,'CA2','CA2-CA3a')
addtetrodelocation_kk('/data12/mari/Egy/','egy',18:21,'ctx','ctx-super-overhpc')
% todo depths

%taskstruct
createtaskstruct('/data12/mari/Egy/','egy',[1 2; 1 4; 1 7],'getcoord_lineartrack');
createtaskstruct('/data12/mari/Egy/','egy',[2 2; 2 4],'getcoord_lineartrack');
createtaskstruct('/data12/mari/Egy/','egy',[3 2; 3 4; 3 6],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[4 2; 4 4; 4 6],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[5 2; 5 4; 5 6],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[6 2; 6 4; 6 7],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[7 2; 7 4; 7 6],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[8 2; 8 4; 8 6],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[9 2; 9 4; 9 6],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[10 2; 10 4; 10 6],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[11 2; 11 4; 11 6],'getcoord_wtrack');
createtaskstruct('/data12/mari/Egy/','egy',[12 2; 12 4; 12 6],'getcoord_wtrack');
kk_updatetaskstruct_sleep('/data12/mari/Egy/','egy',1:12)
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[1 2; 1 4; 1 6],'LinearA')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[2 2; 2 4],'LinearA')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[3 2; 3 4 ; 3 6],'WTrackA')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[4 2; 4 4 ; 4 6],'WTrackA')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[5 2; 5 4 ; 5 6],'WTrackA')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[6 2],'WTrackA')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[6 4 ; 6 7],'WTrackB')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[7 2; 7 4 ; 7 6],'WTrackB')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[8 2; 8 4 ; 8 6],'WTrackB')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[9 2; 9 4 ; 9 6],'WTrackB')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[10 2; 10 4 ; 10 6],'WTrackB')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[11 2; 11 4 ; 11 6],'WTrackB')
kk_updatetaskstruct_environment('/data12/mari/Egy/','egy',[12 2; 12 4 ; 12 6],'WTrackB')
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[1 2 inf; 1 4 inf; 1 7 inf])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[2 2 inf; 2 4 inf])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[3 2 1; 3 4 2; 3 6 3])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[4 2 4; 4 4 5; 4 6 6])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[5 2 7; 5 4 8; 5 6 9])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[6 2 10; 6 4 1; 6 7 2])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[7 2 3; 7 4 4; 7 6 5])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[8 2 6; 8 4 7; 8 6 8])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[9 2 9; 9 4 10; 9 6 11])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[10 2 12; 10 4 13; 10 6 14])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[11 2 15; 11 4 16; 11 6 17])
kk_updatetaskstruct_exposure('/data12/mari/Egy/','egy',[12 2 18; 12 4 19; 12 6 20])

% Dereferencing: eeggnd

addgndtoeeg('egy',1,1:8,1:21)
addgndtoeeg('egy',2,1:7,1:21)
addgndtoeeg('egy',3,1:7,1:21)
addgndtoeeg('egy',4,1:7,1:21)
addgndtoeeg('egy',5,1:7,1:21)
addgndtoeeg('egy',6,1:8,1:21)
addgndtoeeg('egy',7,1:8,1:21)
addgndtoeeg('egy',8,1:8,1:21)
addgndtoeeg('egy',9,7,17)
addgndtoeeg('egy',10,1:8,1:21)
addgndtoeeg('egy',11,1:7,1:21)
addgndtoeeg('egy',12,1:8,1:21)

    



%% DAVE

cd('/data12/kkay/DaveData/')
kk_dayprocess_2('dave01','/data12/kkay/Dav','dav',1,'cmperpix',[.2350 .5 .2350 .5 .2350 .5 .2350])
cd('/data12/kkay/DaveData/')
kk_dayprocess_2('dave02','/data12/kkay/Dav','dav',2,'cmperpix',[.2350 .9 .2350 .9 .2350 .9 .2350])
cd('/data12/kkay/DaveData/')
kk_dayprocess_2('dave03','/data12/kkay/Dav','dav',3,'cmperpix',[.2350 .9 .2350 .9 .2350 .9 .2350])
cd('/data12/kkay/DaveData/')
kk_dayprocess_2('dave04','/data12/kkay/Dav','dav',4,'cmperpix',[.2350 .9 .2350 .9 .2350 .9 .2350])
cd('/data12/kkay/DaveData/')
kk_dayprocess_2('dave05','/data12/kkay/Dav','dav',5,'cmperpix',[.2350 .2350 .9 .2350 .95 .2350])
cd('/data12/kkay/DaveData/')
kk_dayprocess_2('dave06','/data12/kkay/Dav','dav',6,'cmperpix',[.2350 .9 .2350 .9 .2350 .95 .2350])
cd('/data12/kkay/DaveData/')
kk_dayprocess_2('dave07','/data12/kkay/Dav','dav',7,'cmperpix',[.2350 .5 .5 .2350 .5 .2350 .5 .2350 .95])
cd('/data12/kkay/DaveData/')
kk_dayprocess_2('dave08','/data12/kkay/Dav','dav',8,'cmperpix',[.2350 .5 .5 .2350 .5 .2350])








                         

%% TASKSTRUCT

% (after having run taskstruct and lineardayprocess, and
% kk_updatetaskstruct_sleep)


%Dave
createtaskstruct('/data12/kkay/Dav/','dav',[1 2; 1 4; 1 6],'getcoord_lineartrack');
createtaskstruct('/data12/kkay/Dav/','dav',[2 2; 2 4; 2 6],'getcoord_wtrack');
createtaskstruct('/data12/kkay/Dav/','dav',[3 2; 3 4; 3 6],'getcoord_wtrack');
createtaskstruct('/data12/kkay/Dav/','dav',[4 2; 4 4; 4 6],'getcoord_wtrack');
createtaskstruct('/data12/kkay/Dav/','dav',[5 3; 5 5],'getcoord_wtrack');
createtaskstruct('/data12/kkay/Dav/','dav',[6 2; 6 4; 6 6],'getcoord_wtrack');
createtaskstruct('/data12/kkay/Dav/','dav',[7 2; 7 3; 7 5; 7 7; 7 9],'getcoord_wtrack');
createtaskstruct('/data12/kkay/Dav/','dav',[8 2; 8 3; 8 5],'getcoord_wtrack');
kk_updatetaskstruct_sleep('/data12/kkay/Dav/','dav',1:8)
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[1 2; 1 4; 1 6],'LinearA')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[2 2; 2 4; 2 6],'WTrackA')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[3 2; 3 4; 3 6],'WTrackA')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[4 2; 4 4; 4 6],'WTrackA')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[5 3],'WTrackA')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[5 5],'WTrackB')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[6 2; 6 4],'WTrackA')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[6 6],'WTrackB')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[7 2; 7 3; 7 5; 7 7],'WTrackA')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[7 9],'WTrackB')
kk_updatetaskstruct_environment('/data12/kkay/Dav/','dav',[8 2; 8 3; 8 5],'WTrackA')
%cell info, tetinfo
createcellinfostruct('/data12/kkay/Dav/','dav')
createtetinfostruct('/data12/kkay/Dav/','dav')
 % todo locations
 % todo depths
 
 
%% BON

%dereference

addgndtoeeg('bon',3,1:7,[1:8 10:15 17:25 27:30])
addgndtoeeg('bon',4,1:7,[1:8 10:15 17:25 27:30])
addgndtoeeg('bon',5,1:9,[1:8 10:15 17:25 27:30])
addgndtoeeg('bon',6,1:7,[1:8 10:15 17:25 27:30])
addgndtoeeg('bon',7,1:7,[1:8 10:15 17:25 27:30])
addgndtoeeg('bon',8,1:7,[1:8 10:15 17:25 27:30])
addgndtoeeg('bon',9,1:7,[1:8 10:15 17:25 27:30])
addgndtoeeg('bon',10,1:7,[1:8 10:15 17:25 27:30])

%hemisphere field
tettocellinfo('/data12/kkay/Bon/','bon',1:29,'hemisphere')


%% FRA

%dereference

addgndtoeeg('fra',1,[2 4],[1:15 17:25 27:30])
addgndtoeeg('fra',2,[1:7],[1:15 17:25 27:30])
addgndtoeeg('fra',3,[1:7],[1:15 17:25 27:30])  % non-interpolable discrepancy errors noted for tetrode 18 and above, noted 6.10.13 
addgndtoeeg('fra',4,[2:7],[1:15 17:25 27:30])  % non-interpolable discrepancy errors noted for tetrode 18 and above, noted 6.10.13 
addgndtoeeg('fra',5,[1:7],[1:15 17:25 27:30])  % non-interpolable discrepancy errors noted for tetrode 18 and above, noted 6.10.13
addgndtoeeg('fra',6,[1:7],[1:15 17:25 27:30])
addgndtoeeg('fra',7,[1:7],[1:15 17:25 27:30])  % non-interpolable discrepancy errors noted for tetrode 18 and above, noted 6.10.13
addgndtoeeg('fra',8,[1:7],[1:15 17:25 27:30])  % non-interpolable discrepancy errors noted for tetrode 18 and above, noted 6.10.13
addgndtoeeg('fra',9,[1:7],[1:15 17:25 27:30])
addgndtoeeg('fra',10,[1:7],[1:15 17:25 27:30])
addgndtoeeg('fra',11,[1:7],[1:15 17:25 27:30])
addgndtoeeg('fra',12,[1:7],[1:15 17:25 27:30]) % non-interpolable discrepancy errors noted for tetrode 18 and above, noted 6.10.13


% hemisphere field
tettocellinfo('/data12/kkay/Fra/','fra',1:29,'hemisphere')

%% COR

tettocellinfo('/data12/kkay/Bon/','bon',1:29,'hemisphere')
kk_updatetaskstruct_sleep('/data12/kkay/Cor','Cor',1:9)

%dereference
for d = 1:15
    for e=1:15
        for t=1:30
            try
            addgndtoeeg('Cor',d,e,t)
            catch
                continue
            end
        end
    end
end


%% ARN (Annabelle's)


%dereference

addgndtoeeg('arn',5,1:6,1:8)
addgndtoeeg('arn',6,1:7,1:8)
    % bypassed days 11 and 12, odd mismatch in eeg file for day 11 epoch 8 tet 4

%% BAR (Annabelle's)

%dereference
for d = 1:30
    for e=1:30
        for t=1:30
            try
            addgndtoeeg('bar',d,e,t)
            catch
                continue
            end
        end
    end
end


%% CAL (Annabelle's)


%dereference
for d = 1:30
    for e=1:30
        for t=1:30
            try
            addgndtoeeg('cal',d,e,t)
            catch
                continue
            end
        end
    end
end


%% DWI (Annabelle's)


%dereference
for d = 1:30
    for e=1:30
        for t=1:30
            try
            addgndtoeeg('dwi',d,e,t)
            catch
                continue
            end
        end
    end
end

