function runEEG(sepNovel, mtm, selectid)
%function runEEG(sepNovel [, mtm, selectid])
%
% Calculate EEG power spectral density
%  sepNovel     0: calc power everywhere, 
%               1: only in novel arm, (either outer arm if fam conf)
%               2: only in fam arm (either outer arm if fam conf)
%  mtm          whether to use multi-taper method

addpath /home/chengs/theta/common

rerun= 1;
%maxEEGAmp= 10; % max. EEG amplitude
minWinLen= 1; % [sec] min. window length
% definition of outer arm
%pad= 5;
%armlength= 65;

if nargin<2; mtm= 1; end
if nargin<3; selectid= 'CA1PE'; end

setRoot;
load(sprintf('%s/data/tetrodes-%s', root, selectid));
tet= tetrodes;
ntet= length(tet.rat);

switch sepNovel
case 0
    novelStr= '';
case 1
    if strfind(selectid, 'Arm') 
        novelStr= '-novelArm';
    else
        novelStr= '-outerArm';
    end
case 2
    if strfind(selectid, 'Arm') 
        novelStr= '-famArm';
    else
        novelStr= '-outerArm';
    end
otherwise
    error('do not know what to do');
end

starttet=1;
%starttet=110;
%starttet=3;

fprintf(1, '\nAnalyzing EEG for %s, novel selection: %s...\n', selectid, novelStr);
monitorProgress(1, ntet);
for itet=starttet:ntet
    monitorProgress(itet, ntet);
    rat= tet.rat{itet};
    d= tet.num(itet,1); e= tet.num(itet,2); t= tet.num(itet,3);
    fprintf(1, '%s [%d %d %d], itet= %d\n', rat,d,e,t,itet);

    if itet==starttet | ~strcmp(rat, tet.rat{itet-1}) | d~= tet.num(itet-1,1)
        load(sprintf('%s/%s/data2/info',root,rat));
        bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
        load(bfile);
        if(mtm)
            psfile= sprintf('%s/%s/results/EEG/ps%.2d%s.mat',root,rat,d,novelStr);
        else
            psfile= sprintf('%s/%s/results/pwelch/ps%.2d%s.mat',root,rat,d,novelStr);
        end
        if(exist(psfile,'file'))
            load(psfile);
        else
            ps= {};
        end
    end

    if ~rerun;
        if ~isempty(ps) & ...
            length(ps)>=d & ~isempty(ps{d}) & ...
            length(ps{d})>=e & ~isempty(ps{d}{e}) & ...
            length(ps{d}{e})>=t & ~isempty(ps{d}{e}{t}) 
                continue;
        end
    end

    fname= sprintf('%s/%s/data/EEG/%seeg%.2d-%d-%.2d.mat',root,rat,rat,d,e,t);
    load(fname);
    eval(sprintf('eegstruct= %seeg{%d}{%d}{%d};', rat, d, e, t));
    eeg= eegstruct.data;
    fs= eegstruct.samprate;
    mint= eegstruct.starttime;
    maxt= (length(eegstruct.data)-1)/fs+mint;
    sigma= floor(minWinLen*fs/18);
    kernel= gaussian(sigma, 6*sigma);
    cutoff= 3.6559e-04; % corresponding to a time window of ~250ms

    time= behavdata{d}{e}.time;
    dt= mean(diff(time));

    % find valid times for psd analysis
    if ~sepNovel
        ind= find(behavdata{d}{e}.traj>0);
    else
        % only timesteps on valid traj
        if tet.day(itet)<7;  % novel exposure
            % convert between novel arm and traj
            if tet.newarm(itet) < 5;
                if sepNovel==1 
                    trajs= [0, 1]; % novel trajs
                else
                    trajs= [2, 3]; % fam trajs
                end
            elseif tet.newarm(itet) > 5;
                if sepNovel==1 
                    trajs= [2, 3]; % novel trajs
                else
                    trajs= [0, 1]; % fam trajs
                end
            end
        else  % familiar environment
            trajs= [0:3];
        end
        ind= find(ismember(behavdata{d}{e}.traj,trajs));

%        minx= info{d}{e}.centerlinpos+pad
%        maxx= minx+armlength
        %%@@
        minx= info{d}{e}.centerlinpos;
        maxx= info{d}{e}.maxlinpos;


        % only on outer arm
        ind= ind(behavdata{d}{e}.linpos(ind)>=minx & behavdata{d}{e}.linpos(ind)<=maxx);
    end
    [lo,hi]= findcontiguous(ind);

    % only use windows with minimum size
    winall= [time(lo), time(hi)+dt];
    winall(winall<mint)= mint;
    winall(winall>maxt)= maxt;
    win_valid_ind= diff(winall,[],2)>minWinLen;
    winsel= winall(win_valid_ind,:);
    if isempty(winsel); warning('no valid time windows found'); end
    lo= lo(win_valid_ind);
    hi= hi(win_valid_ind);

    % convert to indices for EEG time
    lt= floor((time(lo)-mint)*fs);
    lt(lt<1)= 1;
    ht= floor((time(hi)-mint)*fs);
    ht(ht>length(eeg))= length(eeg);

    % clean-up, 
    large= 10*std(abs(eeg));

    for i=1:size(winsel,1)
        range= lt(i):ht(i);
        ind= find(abs(eeg(range))>large); 
        if ~isempty(ind)
            invalid= zeros(size(range));
            invalid(ind)= 1;
            invalidF=filtfilt(kernel,1,invalid);
            invalid(invalidF>cutoff)= 1;
            ind= find(invalid);
    if(0)
            figure(1); clf
            subplot(2,2,1)
            plot(range/fs+mint, eeg(range));
            hold on
            plot(range(ind)/fs+mint,eeg(range(ind)), 'r');
            xlabel('time (s)'); ylabel('EEG');
            subplot(2,2,2)
            plot(behavdata{d}{e}.time(lo(i):hi(i)), behavdata{d}{e}.linpos(lo(i):hi(i)))
            xlabel('time (s)'); ylabel('lin pos (cm)');
            hold off
%            subplot(2,2,3)
%            plot(range/fs+mint,[invalid', invalidF'])
            subplot(2,2,4)
            plot(behavdata{d}{e}.time(lo(i):hi(i)), behavdata{d}{e}.traj(lo(i):hi(i)))
            xlabel('time (s)'); ylabel('traj');

            fprintf(1,'Checkout: iwin= %d, press key', i)
            pause;
            fprintf(1,'...\n')
    end
            lt(i)= nan;
            ht(i)= nan;
            [lotmp,hitmp]= findcontiguous(range(invalid==0));
            for j=1:length(lotmp);
                lt(end+1)= lotmp(j);
                ht(end+1)= hitmp(j);
            end

%            subplot(2,2,3)
%            hold on
%            for j=1:length(lotmp);
%                plot([lotmp(j):hitmp(j)]/fs+mint,0, 'g.');
%            end
%            pause;
        else
            % plot pretty EEG snippet
            if 1 
                figure(1); clf
                plot((range-range(1))/fs, eeg(range)/max(eeg(range)), 'r');
                xlabel('time (s)'); ylabel('EEG');
                axis tight
                set(gca, 'xlim', [0,1]);
%                myprint('mini', 'EEG-snippet');
fprintf(1, 'key ...');
pause
fprintf(1, '\n');
            end
        end
    end
    lt= sort(lt(isfinite(lt))); ht= sort(ht(isfinite(ht)));

    winsel= [lt ht]/fs+mint;
    win_valid_ind= diff(winsel,[],2)>minWinLen;
    winsel= winsel(win_valid_ind,:);

    if isempty(winsel); warning('no valid time windows found'); end

    % normalize
    EEG_valid= makecontiguous(lt,ht);
    eegstruct.data= eeg/ sqrt(sum(eeg(EEG_valid).^2));

%    plot(eegstruct.data(EEG_valid));
%    keyboard

    if(mtm)
        [w, f, px] = swmts(eegstruct, winsel);
        ps{d}{e}{t}.winall= winall;
        pxs= mean(px,2);
        ps{d}{e}{t}.win= w;
        ps{d}{e}{t}.freq= f;
        ps{d}{e}{t}.psd= px;
    else % use Welch method to compute power spectrum
        fs= eegstruct.samprate;
        [px,f]=pwelch(eeg,[],[],[],fs);

        ps{d}{e}{t}.freq= f;
        ps{d}{e}{t}.psd= px;
    end
    
    save(psfile, 'ps');
end

