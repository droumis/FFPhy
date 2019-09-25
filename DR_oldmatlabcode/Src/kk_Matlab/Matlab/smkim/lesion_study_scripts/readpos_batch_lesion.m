
% read all pos files and extract them into matlab data files

%subjects = { 'M02' 'M03' 'M06' };
subjects = { 'M12' 'M13' 'M14' 'M16' 'M17' 'M19' 'M20' 'M22' 'M24' 'M25' 'M26' };
days = 1:2;
sessions = 1:2;

%{
load('sessions.mat');
for i = 1:length(sessions)
    for j = 1:length(sessions{i})

    session = sessions{i}{j}.descript;
    environment = sessions{i}{j}.environment;
    try
        posfilename = ['day' sprintf('%1d',i) '_' session '.pos'];
        rawpos{i}{j} = readpos(posfilename, ...
            sessions{i}{j}.tstart,sessions{i}{j}.tend);
        rawpos{i}{j}.subject = subject;
        rawpos{i}{j}.day = i;
        rawpos{i}{j}.session = session;
        rawpos{i}{j}.environment = environment;
    catch
        warning('unable to read %s',posfilename);
    end
    try
        pngfilename = ['day' sprintf('%1d',i) '_' environment '.png'];
        sampleframe = imread(pngfilename);
        rawpos{i}{j}.sampleframe = sampleframe;
    catch
        warning('unaable to read %s',pngfilename);
    end
end
rawposfilename = 'rawpos.mat';
save(rawposfilename,'rawpos');
clear('rawpos');
%}

fid = fopen('round3_times.txt');
for s = 1:length(subjects)
    for i = 1:2
        for j = 1:2
            if i == 1
                posfilename = [subjects{s} '_linear_preop_' num2str(j) '.pos'];
            elseif i == 2
                posfilename = [subjects{s} '_linear_postop_' num2str(j) '.pos'];
            end
            frewind(fid);
            while 1
                nextline = fgetl(fid);
                if nextline == -1
                    break
                elseif regexp(nextline,['^' posfilename])
                    splits = regexp(nextline,'\s*','split');
                    tstart = splits{2}; tend = splits{3};
                end
            end
            try
                rawpos{i}{j} = readpos_lesion(posfilename,tstart,tend);
                rawpos{i}{j}.subject = subjects{s};
                rawpos{i}{j}.environment = 'linear track';
                if i == 1
                    rawpos{i}{j}.day = 'pre surgery';
                elseif i == 2
                    rawpos{i}{j}.day = 'post surgery';
                end
                if j == 1
                    rawpos{i}{j}.session = 'run1';
                elseif j == 2
                    rawpos{i}{j}.session = 'run2';
                end
                if i == 1
                    pngfilename = 'round3_linear_preop.png';
                elseif i == 2
                    pngfilename = 'round3_linear_postop.png';
                end
                sampleframe = imread(pngfilename);
                rawpos{i}{j}.sampleframe = sampleframe;
            end
        end
    end
    rawposfilename = [subjects{s} '_linear_rawpos.mat'];
    try
        save(rawposfilename,'rawpos');
        disp(subjects{s});
    end
    clear('rawpos');
end

