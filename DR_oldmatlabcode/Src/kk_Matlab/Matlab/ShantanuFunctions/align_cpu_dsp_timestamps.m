function postimestamp_align = align_cpu_dsp_timestamps(cpudsp_filename,cpupos_filename,save_filename)

%This function should be used when there is an error in the possynctimes.

%It reads in the cpudsptimecheck and the cpupostimestamp files, looks
%for differences in the clocks and creates a postimestamp file that
%accounts for drift. Note, this is only a good option if the normal
%nspike_postimestamp approach fails.

%Inputs
%cpudsp_filename:   Name of cpudsptimecheck file to use
%cpupos_filename:   Name of cpupostimestamp file to use
%save_filename:     Filename if wish to save output, Default is not to save

if nargin == 2
    save = 0;
else
    save = 1;
end

%Load cpudsp timestamps
fid = fopen(cpudsp_filename);
tf = 0;
while tf == 0
    s = fgets(fid);
    tf = strncmp(s,'%%ENDHEADER',5);
end

cpudsp = fread(fid,'uint32');

%Load cpupos timestamps
fid = fopen(cpupos_filename);
tf = 0;
while tf == 0
    s = fgets(fid);
    tf = strncmp(s,'%%ENDHEADER',5);
end

cpupos = fread(fid,'uint32');

%Determine offset between cpu and dsp clocks
offset = cpudsp(1:3:end) - cpudsp(2:3:end);

%Determine which offset to use for each cpupostimestamp
index = lookup(cpupos,cpudsp(1:3:end));
postimestamp_align = cpupos - offset(index);

if save == 1
    fid = fopen(save_filename,'a+');
    fprintf(fid,'%%%%BEGINHEADER\n');
    fprintf(fid,'%% File type:	Binary\n');
    fprintf(fid,'%% Extraction type:	dsp position frame time stamps\n');
    fprintf(fid,'%% Fields:	  timestamp (unsigned int)\n');
    fprintf(fid,'%%%%ENDHEADER\n');
    
    fwrite(fid, postimestamp_align, 'uint32');
    fclose(fid);
end

end