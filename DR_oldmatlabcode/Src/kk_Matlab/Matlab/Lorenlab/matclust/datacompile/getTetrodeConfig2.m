function output = getTetrodeConfig()

files = dir('*.dat.config');
tetnum = [];
thresh = [];
gain = [];
for i = 1:length(files)
    fid = fopen(files(i).name,'r');
    text = fread(fid,'char');
    text = char(text)';
    sysloc = strfind(text,'datatype')+10;
    spikemachine = strcmp(text(sysloc:sysloc+4),'SPIKE');
    
    if (spikemachine)
        
        tetnumloc = (strfind(text,'number')+7);
        %tetnumloc = tetnumloc(2:2:end);
        %gainloc = (strfind(text,'gain')+5);
        %filterloc = (strfind(text,'filter')-1);
        dispvalloc = (strfind(text,'maxdispval')+11);
        threshloc = (strfind(text,'thresh')+7);
        threshloc = threshloc(2:end);
        tmptetnum = [];
        tmpthresh = [];
        %tmpgain = [];
        for j = 1:length(tetnumloc)
            
            tmptetnum(j) = str2num(text(tetnumloc(j):tetnumloc(j)+1));          
            %j
            %threshloc
            %dispvalloc
            %text(threshloc(j):dispvalloc(j)-15)
            tmpthresh(j) = str2num(text(threshloc(j):dispvalloc(j)-15));
            %tmpgain(j) =  str2num(text(gainloc(j):filterloc(j)));   
        end
        tetnum = [tetnum; tmptetnum'];
        thresh = [thresh; tmpthresh'];
        %gain = [gain; tmpgain'];
        %tetnum =  str2num(text(strfind(text,'number')+7:strfind(text,'number')+8)')
        %tetnum = tetnum(2:2:end)
    end
    
    fclose(fid);
end    

infomatrix = [tetnum thresh];
infomatrix = sortrows(infomatrix,1);

count = 1;
while (count<size(infomatrix,1))
    
    %infostruct(infomatrix(count,1)).gain = infomatrix(count:count+3,3)';
    infostruct(infomatrix(count,1)).thresh = infomatrix(count:count+3,2)';
    count = count+4;
end

output = infostruct;
    
