function selectCells_Adaptx 
  
epoch= 2; 

load /bach/Ter/data/Adaptx3.0t0.05d/teradapt.mat
valcn= find(teradapt(epoch).cellnum(:,1)~=-1);
ci=valcn(1:8:size(valcn,1));
cellselect= teradapt(epoch).cellnum(ci,:);

save /bach/Ter/data2/cellselect-Adaptx3-0t0-05d cellselect