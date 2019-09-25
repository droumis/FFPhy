function []=Raw2PLX(directoryname,Block,animalname,channels,downsample)

% converts TDT data to Plexon binary file
% Block      - data block to be converted eg 20
% animalname - name of animal for processing eg CML9
% channels   - channels to be converted eg [1:64] [1 3 34]
% downsample - take every nth datapoint eg 10


% sampling rate of TDT data is 24414.1Hz
% time step is 0.0000409599370855366s


% modified by JY from Jiwei_Raw2PLX.m

% Block to retrieve
% Block = 44;
iDate = '020110712';
cd(directoryname);
dMaxV = 5e-3; %from Serge's original file - for scaling purposes
iResolutionBits = 16;
dVoltsToIntScale = pow2(iResolutionBits-1)/dMaxV;

chanlist=channels';

timelist=load(sprintf('/data14/jai/%s/timelist.mat',animalname));

timelist=timelist.t;



title(sprintf('CML9 Block %d',Block));

for i=1:size(chanlist,1);
    

    for chan = chanlist(i,1);
        %Open file for writing
        fn = sprintf('/data14/jai/%s/Block-%d/Jai_Block-%d_xWav_ch%d.sev',...
            animalname,Block,Block,chan);
        fmt = 'float32';
        fid = fopen(fn, 'rb');
    %     tmp = fread(fid, 1e9, fmt);
    %     if chan == 1
    %         y = zeros(length(tmp),length(chan));
    %     end
    %     y (:,chan)= tmp 

        y = fread(fid, 1e9, fmt);
        y = y * dVoltsToIntScale;

        mat=[];
        mat = y(1:downsample:end,1);
        timestepms = 1000*0.0000409599370855366*(downsample);

        fclose(fid);

        %where data is stored

        direct=sprintf('/data14/jai/%s/',...
            animalname);
        cd(direct);
        datadir = dir('Plexon');

            if (isempty(datadir));
                %an a plot folder needs to be created
                !mkdir Plexon
            end

        fsave = sprintf ('/data14/jai/%s/Plexon/Block-%d_%d_C%d-%s.mat', ...
            animalname,Block, chan, iDate);
        % disp(['read ' num2str(size(y,1)) ' elements']);
        % f = fopen(fsave, 'w+b');
        %save(fsave,'mat','timestepms'); % convert to int16 and write
        % fclose(f);



        timefile=sprintf('/data14/jai/%s/LaserStim_Params_Jai_Blk%d.mat',...
            animalname,Block);
        load(timefile);

        timestep=[0:timestepms/1000:timestepms*size(mat,1)/1000]';
        
        %[filtwts] = designeegfilt(1/(0.0000409599370855366*downsample),300,6000)
        
       
        %mat = filtfilt(filtwts,1,mat);
        
        mat=[mat timestep(1:end-1,1)];
        
        mat=mat(1:100:end,:);
        
        figure;
   
        subplot(size(chanlist,1),1,i);
        
        title(sprintf('10x 1s 594nm'));
        
        for u=1:size(timelist,1);
            subplot(size(timelist,1),1,u);
            figurename=sprintf('CML9_B_%d_C_%d.pdf',Block,channels);
            
            %set(gca,'ytick',);
            yellow=[255 192 0]/255;
            
            stim(:,1)=Data.Pulse(:,1);
            stim(:,2)=Data.Pulse(:,1)+Data.PulseDur(:,1)/1000-0.05;
            
            for i =1:size(stim,1);
                line([stim(i,1) stim(i,2)],[-2000 -2000],'Color',yellow,'LineWidth',5000);
              
            end
            
            hold on;
          
            plot(mat(:,2),mat(:,1));
            
            hold on;

            
           
            ylabel(sprintf('%sV',num2str(timelist(u,1),1)),'Rotation',0);
            ylim([-10e3 10e3]);
            
            xl=[timelist(u,3)-2 timelist(u,3)+12];
            set(gca,'xlim',xl);
            
            
            
            u=u+1;
        end
        
        set(gca,'xtick',[],'ytick',[]);
        filename=sprintf('%s_%s_%s_LFP.pdf',animalname,num2str(Block),num2str(chan));
        saveas(gcf,filename,'pdf');
        
        close;

    
end

set(gca,'xtickmode','auto');

end

% clear f y fid fsave dMaxV iResolutionBits ;
% clear iDate fn fmt tmp dVoltsToIntScale;
