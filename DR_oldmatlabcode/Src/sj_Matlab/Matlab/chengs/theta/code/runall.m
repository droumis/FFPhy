function runall(indir, inSelId, eps_spatial, eps_temp, rerun)
%function runall
%function runall(indir)
%function runall(indir, inSelId)
%function runall(indir, inSelId)
%function runall(indir, inSelId, range_spatial, range_temp)

global REDO rats

%rats={'kyl', 'ter', 'sta', 'fel'};
rats={'gen'};

if nargin<5; REDO= 0; else REDO= rerun; end


if nargin>=1  | isempty(indir)
    if ~iscell(indir)
        adapt_dir= {indir};
    else
        adapt_dir= indir;
    end
else
    adapt_dir= { 'nonripple_xp3_t01'; };
end

if nargin<2 | isempty(inSelId)
    select_id= { 'pf9-novel', 'pf9-fam',
        %'pf9-novelArm', %'pf9-famArm',
    };
else
    if ~iscell(inSelId)
        select_id= {inSelId};
    else
        select_id= inSelId;
    end
end

if nargin>=4
    if iscell(indir); error('can only handle one directory prefix'); end
    adapt_dir= {};
    nx= length(eps_spatial);
    nt= length(eps_temp);
    disp('Generating directory names...');
    for ix=1:nx
        eps_x= eps_spatial(ix);
        for it=1:nt
            eps_t= eps_temp(it);
%            dirname= auxGetResultsDir(indir, eps_x, eps_t);
            dirname= auxMakeResultsDir(indir, eps_x, eps_t);
            adapt_dir{end+1}= dirname;
        end
    end
end

for i=1:length(adapt_dir)
    for iS=1:length(select_id)
        fprintf(1, '******* starting %s %s %s *******\n', ...
            select_id{iS}, adapt_dir{i}, datestr(clock));
        t0= clock;
        switch REDO
        case {-1}
        case {0} % no reruns
            forallrats('runAdaptFilter([], [], 1)', select_id{iS}, adapt_dir{i}, rats);
        case {1}     % run all 
            forallrats('runAdaptFilter', select_id{iS}, adapt_dir{i}, rats);
        end

        forallrats('showKSplot([],1,r);', select_id{iS}, adapt_dir{i}, rats);
        fprintf(1, 'KS test %s %s: ', select_id{iS}, adapt_dir{i});
        showKSplot([],1,-1); % show summary across all animals

%         get stats, using common stat options
%        forallrats('runStats2', select_id{iS}, adapt_dir{i}, rats);

        fprintf(1, '******* done with %s %s %s *******\n', ...
            select_id{iS}, adapt_dir{i}, datestr(clock));
        et= etime(clock, t0);
        eh= floor(et/3600); em= floor((et-3600*eh)/60);  es= et-3600*eh-60*em;
        fprintf(1, 'runtime: %d h %d min %.1f sec\n', eh, em, es);
    end
end

function str= auxGetFilenameString(x);

if x== 0
    str= '0';
elseif x>= 0.001
    str= num2str(x);
    if x < 1
        i= strfind(str, '.');
        if ~isempty(i)
            str= [str(1:i-1) str(i+1:end)];
        end
    end
else
    str= num2str(x, '%.0e');
    i= strfind(str, 'e-0');
    str= [str(1:i+1) str(i+3:end)];
end

function dirname= auxGetResultsDir(indir, eps_x, eps_t)

xstr= auxGetFilenameString(eps_x);
tstr= auxGetFilenameString(eps_t);

dirname= sprintf('%s%s_t%s', indir, xstr, tstr);


function dirname= auxMakeResultsDir(indir, eps_x, eps_t)

global REDO rats

xstr= auxGetFilenameString(eps_x);
tstr= auxGetFilenameString(eps_t);

dirname= sprintf('%s%s_t%s', indir, xstr, tstr);

root= '/home/chengs/theta';

for i=1:length(rats)
    name= [root '/' rats{i} '/results/' dirname]; 
    fprintf(1, '%s ',  name);
    if exist(name)==7 & REDO~= 1
        fprintf(1, ' -- using existing files.\n');
    else
        if exist(name)==7 
            fprintf(1, ' -- overwriting files.\n');
        else
            fprintf(1, ' -- generating new files.\n');
        end

        if eps_x==0
            fid= auxCopyDir(rats{i}, dirname, 'getLocalFilterModelConst_ISI.m');
            fprintf(fid, 'fopts.eps= [ %g * ones(1,4)];\n', eps_t);


        else
            switch indir(end-1:end)
            case '_x'
                if eps_t > 0
                    fid= auxCopyDir(rats{i}, dirname, 'getLocalFilterModel1d.m');
                    fprintf(fid, 'fopts.eps= [ %g  * ones(1,4), %g * ones(1,4)];\n\n', eps_x, eps_t);
                else
                    fid= auxCopyDir(rats{i}, dirname, 'getLocalFilterModel1d_Const.m');
                    fprintf(fid, 'fopts.eps= [ %g * ones(1,4) ];\n\n', eps_x);
                end
            case 'xp'
                if eps_t > 0
                    fid= auxCopyDir(rats{i}, dirname, 'getLocalFilterModel2d.m');
                    fprintf(fid, 'fopts.eps= [ %g  * ones(1,16), %g * ones(1,4)];\n\n', eps_x, eps_t); 
                else
                    fid= auxCopyDir(rats{i}, dirname, 'getLocalFilterModel2d_Const.m');
                    fprintf(fid, 'fopts.eps= [ %g  * ones(1,16) ];\n\n', eps_x);
                end
            otherwise 
                error('do not know what to do');
            end
        end
        if eps_x > 0
            if eps_x <= 0.05
                fprintf(fid, 'model.spatial.conv= %g;\n', eps_x/10);
            else
                fprintf(fid, 'model.spatial.conv= %g;\n', .005);
            end
        end

        if eps_t > 0
            if eps_t <= 0.05
                fprintf(fid, 'model.isi.conv= %g;\n', eps_t/10);
            else
                fprintf(fid, 'model.isi.conv= %g;\n', .005);
            end
        end

        if strfind(name, 'longISI')
            fprintf(fid, '\nmodel.max_isi= 0.5;           %% in sec!\n');
            fprintf(fid, 'model.isi.cpx=[-small 1 3 5 7 11 15:10:485 500+small 501]/1000;\n');
        end

        fclose(fid);
    end
end

function fid= auxCopyDir(rat, dirname, optfile)
global REDO
root= '/home/chengs/theta';
name= [root '/' rat '/results/' dirname]; 
copyfile([root '/' rat '/results/templ'], name);
if REDO==1
    delete([name '/adaptest*.mat']);
end
copyfile([root '/' optfile], [name '/getLocalFilterModel.m'])
fid= fopen([name '/getLocalFilterModel.m'], 'a');

