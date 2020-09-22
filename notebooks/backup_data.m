%% "newer" mode (-u or --update flag) plus -t (to copy file modified time), -r (for recursive folders), and -v (for verbose output to see what it is doing)

%         !rsync -utrvn --exclude '.Trash*' --exclude 'lost+found' -e 'ssh -p 7777' droumis@169.230.191.64:/stelmo/demetris/analysis/  /media/droumis/DR_swapdata11/analysis/ &&
        
%     if filterframework:
%         print('backing up vortex2ToSwap filterframework')
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/D10/filterframework/ /media/droumis/DR_swapdata10/D10/filterframework/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/D12/filterframework/ /media/droumis/DR_swapdata10/D12/filterframework/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/JZ2/filterframework/ /media/droumis/DR_swapdata10/JZ2/filterframework/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/D13/filterframework/ /media/droumis/DR_swapdata10/D13/filterframework/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/JZ3/filterframework/ /media/droumis/DR_swapdata10/JZ3/filterframework/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/JZ1/filterframework/ /media/droumis/DR_swapdata10/JZ1/filterframework/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/JZ4/filterframework/ /media/droumis/DR_swapdata10/JZ4/filterframework/
% 
%     if mountainlab_output:
%         print('backing up vortex2ToSwap mountainlab_output')
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/D10/mountainlab_output/ /media/droumis/DR_swapdata10/D10/mountainlab_output/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/D12/mountainlab_output/ /media/droumis/DR_swapdata10/D12/mountainlab_output/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/JZ2/mountainlab_output/ /media/droumis/DR_swapdata10/JZ2/mountainlab_output/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/D13/mountainlab_output/ /media/droumis/DR_swapdata10/D13/mountainlab_output/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/JZ3/mountainlab_output/ /media/droumis/DR_swapdata10/JZ3/mountainlab_output/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/JZ1/mountainlab_output/ /media/droumis/DR_swapdata10/JZ1/mountainlab_output/
%         !rsync -utrv --exclude '.Trash*' --exclude 'lost+found' /data2/demetris/JZ4/mountainlab_output/ /media/droumis/DR_swapdata10/JZ4/mountainlab_output/


     