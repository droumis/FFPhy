function [ output_args ] = animalinformation( name,strain,sex,DOB,wrstart,trainstart,implant, retrain,recording,lesion,sacrifice,comments)
%Information about rat all dates are YYYMMDD
% name     - name of animal eg. 'H2'
% strain    - strain of rat
% sex       - sex
% DOB       - date of birth
% wrstart   - start date of weight restriction
% trainstart- start date of initial training
% implant   - date of implant surgery
% retrain   - date of retraining after surgery
% recording - date of recording
% lesion    - date of electrolytic lesioning 
% sacrifice - date of sacrifice

animalinformation.name=name;     
animalinformation.strain=strain;    
animalinformation.sex=sex;       
animalinformation.DOB=DOB;      
animalinformation.wrstart=wrstart;   
animalinformation.trainstart=trainstart;
animalinformation.implant = implant; 
animalinformation.retrain = retrain;
animalinformation.recording =recording;
animalinformation.lesion   = lesion;
animalinformation.sacrifice =sacrifice;
animalinformation.comments =comments;

save('animalinformation.mat','animalinformation');



end

