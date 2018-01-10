



%backup rsync data from local SSDs to HDDs on data20 and vortex1
swapToData20 = 0;
swapToVortex1 = 0;
data20ToVortex1 = 0;
vortex1ToVortex = 0;

preprocessing = 0;
filterframework = 0;
analysis = 0;
raw = 0;
%% to data 20
if swapToData20
    if preprocessing
        % !rsync -avp /opt/DR_swapdata4/D10/preprocessing/ /opt/data20/D10/preprocessing/
        % !rsync -avp /opt/DR_swapdata5/D12/preprocessing/ /opt/data20/D12/preprocessing/
        % !rsync -avp /opt/DR_swapdata6/D13/preprocessing/ /opt/data20/D13/preprocessing/
        % !rsync -avp /opt/DR_swapdata8/JZ2/preprocessing/ /opt/data20/JZ2/preprocessing/
        % !rsync -avp /opt/DR_swapdata7/JZ3/preprocessing/ /opt/data20/JZ3/preprocessing/
        % !rsync -avp /opt/DR_swapdata9/JZ4/preprocessing/ /opt/data20/JZ4/preprocessing/
        
        !rsync -avp /opt/DR_swapdata4/ /opt/data20/
        !rsync -avp /opt/DR_swapdata5/ /opt/data20/
        !rsync -avp /opt/DR_swapdata6/ /opt/data20/
        !rsync -avp /opt/DR_swapdata8/ /opt/data20/
        !rsync -avp /opt/DR_swapdata7/ /opt/data20/
        !rsync -avp /opt/DR_swapdata9/ /opt/data20/
    end
    if filterframework
        !rsync -avp /opt/DR_swapdata10/ /opt/data20/
    end
    if analysis
        !rsync -avp /opt/DR_swapdata11/analysis/ /opt/data20/analysis/
    end
end
%% swap to vortex1
if swapToVortex1
    if preprocessing
        % !rsync -avp /opt/DR_swapdata4/D10/preprocessing/ /opt/data20/D10/preprocessing/
        % !rsync -avp /opt/DR_swapdata5/D12/preprocessing/ /opt/data20/D12/preprocessing/
        % !rsync -avp /opt/DR_swapdata6/D13/preprocessing/ /opt/data20/D13/preprocessing/
        % !rsync -avp /opt/DR_swapdata8/JZ2/preprocessing/ /opt/data20/JZ2/preprocessing/
        % !rsync -avp /opt/DR_swapdata7/JZ3/preprocessing/ /opt/data20/JZ3/preprocessing/
        % !rsync -avp /opt/DR_swapdata9/JZ4/preprocessing/ /opt/data20/JZ4/preprocessing/
        
        !rsync -avp /opt/DR_swapdata4/ /opt/vortex1/demetris/
        !rsync -avp /opt/DR_swapdata5/ /opt/vortex1/demetris/
        !rsync -avp /opt/DR_swapdata6/ /opt/vortex1/demetris/
        !rsync -avp /opt/DR_swapdata8/ /opt/vortex1/demetris/
        !rsync -avp /opt/DR_swapdata7/ /opt/vortex1/demetris/
        !rsync -avp /opt/DR_swapdata9/ /opt/vortex1/demetris/
    end
    if filterframework
        !rsync -avp /opt/DR_swapdata10/ /opt/vortex1/demetris/
    end
    if analysis
        !rsync -avp /opt/DR_swapdata11/analysis/ /opt/vortex1/demetris/analysis/
    end
end
%% to vortex1
if data20ToVortex1
    if preprocessing
        %preprocessing
%         !rsync -avp /opt/data20/D10/preprocessing/ /opt/vortex1/demetris/D10/preprocessing/
%         !rsync -avp /opt/data20/D12/preprocessing/ /opt/vortex1/demetris/D12/preprocessing/
%         !rsync -avp /opt/data20/D13/preprocessing/ /opt/vortex1/demetris/D13/preprocessing/
%         !rsync -avp /opt/data20/JZ2/preprocessing/ /opt/vortex1/demetris/JZ2/preprocessing/
%         !rsync -avp /opt/data20/JZ3/preprocessing/ /opt/vortex1/demetris/JZ3/preprocessing/
%         !rsync -avp /opt/data20/JZ4/preprocessing/ /opt/vortex1/demetris/JZ4/preprocessing/
%         !rsync -avp /opt/data20/JZ1/preprocessing/ /opt/vortex1/demetris/JZ1/preprocessing/
        
        !rsync -avp /opt/data20/ /opt/vortex1/demetris/
    end
    if filterframework
%         !rsync -avp /opt/DR_swapdata10/ /opt/vortex1/demetris/
        
        !rsync -avp /opt/data20/D10/filterframework/ /opt/vortex1/demetris/D10/filterframework/
        !rsync -avp /opt/data20/D12/filterframework/ /opt/vortex1/demetris/D12/filterframework/
        !rsync -avp /opt/data20/D13/filterframework/ /opt/vortex1/demetris/D13/filterframework/
        !rsync -avp /opt/data20/JZ2/filterframework/ /opt/vortex1/demetris/JZ2/filterframework/
        !rsync -avp /opt/data20/JZ3/filterframework/ /opt/vortex1/demetris/JZ3/filterframework/
        !rsync -avp /opt/data20/JZ4/filterframework/ /opt/vortex1/demetris/JZ4/filterframework/
        !rsync -avp /opt/data20/JZ1/filterframework/ /opt/vortex1/demetris/JZ1/filterframework/
    end
    if analysis
%         !rsync -avp /opt/DR_swapdata11/analysis/ /opt/vortex1/demetris/analysis/
        !rsync -avp /opt/data20/analysis/ /opt/vortex1/demetris/analysis/
    end
end
%% to vortex0
if vortex1ToVortex
    if raw
        !rsync -avp /opt/vortex1/demetris/D10/raw/ /opt/vortex/demetris/D10/raw/
        !rsync -avp /opt/vortex1/demetris/D12/raw/ /opt/vortex/demetris/D12/raw/
        !rsync -avp /opt/vortex1/demetris/D13/raw/ /opt/vortex/demetris/D13/raw/
        !rsync -avp /opt/vortex1/demetris/JZ1/raw/ /opt/vortex/demetris/JZ1/raw/
        !rsync -avp /opt/vortex1/demetris/JZ2/raw/ /opt/vortex/demetris/JZ2/raw/
        !rsync -avp /opt/vortex1/demetris/JZ3/raw/ /opt/vortex/demetris/JZ3/raw/
    end
end
% !rsync -avp /opt/vortex1/demetris/JZ4/raw/ /opt/vortex/demetris/JZ4/raw/
% %no room on vortex for JZ4 raw.. backup is on rdx

