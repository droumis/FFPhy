%------------------------------------%
%| Save raw as .i16 file for Plexon  |
%------------------------------------%

% Block to retrieve
Block = 44;
iDate = '07-08-2011';

dMaxV = 5e-3; %from Serge's original file - for scaling purposes
iResolutionBits = 16;
dVoltsToIntScale = pow2(iResolutionBits-1)/dMaxV;

for chan = 1:64
    %Open file for writing
    fn = sprintf('E:/RS4/Block-%d/JiWei_Block-%d_xWav_ch%d.sev',...
        Block,Block,chan);
    fmt = 'float32';
    fid = fopen(fn, 'rb');
%     tmp = fread(fid, 1e9, fmt);
%     if chan == 1
%         y = zeros(length(tmp),length(chan));
%     end
%     y (:,chan)= tmp 

    y = fread(fid, 1e9, fmt);
    y = y * dVoltsToIntScale;
    fclose(fid);
    
    %where data is stored
    fsave = sprintf ('E:/Plexon Files/JiWei_%d_C%d-%s.i16', ...
        Block, chan, iDate);
    disp(['read ' num2str(size(y,1)) ' elements']);
    f = fopen(fsave, 'w+b');
    fwrite(f, y, 'int16'); % convert to int16 and write
    fclose(f)
end

clear f y fid fsave dMaxV iResolutionBits ;
clear iDate fn fmt tmp dVoltsToIntScale;
