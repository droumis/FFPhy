function dout= selectData(data,T)

dout.timearray= data.timearray(1:T);
dout.posarray= data.posarray(1:T);
dout.phasearray= data.phasearray(1:T);
dout.fieldID= data.fieldID(1:T);

tmax= data.timearray(T);
sel= find(data.spiketimes < tmax);
dout.spiketimes= data.spiketimes(sel);
dout.ispikes= data.ispikes(sel);

