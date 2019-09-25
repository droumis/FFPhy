function allEEGs=plotEEGs3(pathString,day,epoch,ref)
 if ref==0  
           load([pathString '0' num2str(1)]);
refTrace=eeg{day}{epoch}{1}.data*0;
     

 else
      load([pathString '0' num2str(ref)]);

 refTrace=eeg{day}{epoch}{ref}.data;
     
 end

figure
allEEGs=[];
for tets1=1:21,
    if floor(tets1/10)==0
        load([pathString '0' num2str(tets1)]);
    else
        load([pathString num2str(tets1)]);
    end
    eegTrace=eeg{day}{epoch}{tets1}.data;
        length(eegTrace)

    eegTrace=eegTrace(1:min(length(refTrace),length(eegTrace)));
    if tets1~=ref
           eegTrace=eegTrace(1:min(length(refTrace),length(eegTrace))) + refTrace(1:min(length(refTrace),length(eegTrace)));
    end
    
    allEEGs=[allEEGs [eegTrace; zeros(length(allEEGs)-length(eegTrace),1)]];
%    allEEGs=[allEEGs eegTrace(1:length(allEEGs(:,1)))];
        
    
    if tets1<8
    plot(eegTrace+1500*tets1);hold on;
    elseif tets1<15
    plot(eegTrace+1500*tets1,'r');hold on;
    else
        plot(eegTrace+1500*tets1,'k');hold on;
    
    end
end
