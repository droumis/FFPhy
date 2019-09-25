cellselect.cellnum=[1 2 1 1]
cellselect.traj{1}=[0]
save cellselect-test cellselect

timeselect{1}{2}.time{1}=[100:100:800]
timeselect{1}{2}.index{1}=[1e5:1e5:8e5]
save timeselect-100s timeselect

timeselect{1}{2}.time{1}=[200:200:800]
timeselect{1}{2}.index{1}=[2e5:2e5:8e5]
save timeselect-200s timeselect

timeselect{1}{2}.time{1}=[400:400:800]
timeselect{1}{2}.index{1}=[4e5:4e5:8e5]
save timeselect-400s timeselect

timeselect{1}{2}.time{1}=[800]
timeselect{1}{2}.index{1}=[8e5]
save timeselect-end timeselect
