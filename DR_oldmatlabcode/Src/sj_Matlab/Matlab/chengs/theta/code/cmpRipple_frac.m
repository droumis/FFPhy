N=[];

d=1;
load /bach/theta/work/ripple/data_CA1PE_day1.mat
for ia=1:2; N(d,ia,1)= sum(sum(nRips{ia}{2})); end

d=2;
load /bach/theta/work/ripple/data_CA1PE_day2.mat
for ia=1:2; N(d,ia,1)= sum(sum(nRips{ia}{2})); end

d=3;
load /bach/theta/work/ripple/data_CA1PE_day3.mat
for ia=1:2; N(d,ia,1)= sum(sum(nRips{ia}{2})); end


d=1;
load /bach/theta/work/tripple/data_CA1PE_day1.mat
for ia=1:2; N(d,ia,2)= sum(sum(nRips{ia}{2})); end

d=2;
load /bach/theta/work/tripple/data_CA1PE_day2.mat
for ia=1:2; N(d,ia,2)= sum(sum(nRips{ia}{2})); end

d=3;
load /bach/theta/work/tripple/data_CA1PE_day3.mat
for ia=1:2; N(d,ia,2)= sum(sum(nRips{ia}{2})); end

fprintf(1,'fraction of ripples in quiescence\n    famArm    novelArm\n');

disp(1-N(:,:,2)./N(:,:,1))
