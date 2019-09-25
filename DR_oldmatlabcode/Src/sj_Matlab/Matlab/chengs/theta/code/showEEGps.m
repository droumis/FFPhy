function showEEGps(sepNovel, mtm, arg3, arg4)
%function showEEGps(sepNovel, mtm)
%function showEEGps(sepNovel, mtm, selectid)
%function showEEGps(sepNovel, mtm, rat, num)
%
% Analyze EEG power spectrum.
%  sepNovel     0: calc power everywhere, 
%               1: only in novel arm, (either outer arm if fam conf)
%               2: only in fam arm (either outer arm if fam conf)
%  mtm          whether to use multi-taper method


showPS= 1;
pausePlot= 0;
printPS= 0;

cutoff= 180; % (sec)
%cutoff= -1;


PSfac= 1e7;

setRoot;
if nargin== 4
    tet.rat= {arg3};
    tet.num= arg4;
    ntet= 1;
    selectid= 'pf9-novelArm';
else
    if nargin==3
        selectid= arg3;
    else 
        selectid= 'CA1PE';
    end
    load(sprintf('%s/data/tetrodes-%s', root, selectid));
    tet= tetrodes;
    ntet= length(tet.rat);
%    load(sprintf('%s/data/epochs-%s', root, selectid));
%    ep= epochs;
%    nep= length(ep.rat);
end

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


ana.theta= zeros(1,ntet);
ana.max= zeros(1,ntet);
ana.width= zeros(1,ntet);
if(mtm)
    ana.stdmax= zeros(1,ntet);
    figstr= 'pmtm';
else
    figstr= 'pwelch';
end

starttet= 1;
%starttet= 110;
%starttet= 223;

global lindistpos  

for itet=starttet:ntet
    rat= tet.rat{itet};
    d= tet.num(itet,1); e= tet.num(itet,2); t= tet.num(itet,3);

    % change directory, set local options
    if itet==starttet | ~strcmp(rat, tet.rat{itet-1}) | d~= tet.num(itet-1,1)
        if(mtm)
            psfile= sprintf('%s/%s/results/EEG/ps%.2d%s.mat',root,rat,d,novelStr);
        else
            psfile= sprintf('%s/%s/results/pwelch/ps%.2d%s.mat',root,rat,d,novelStr);
        end
        fprintf(1, 'loading %s ...\n', psfile);
        load(psfile);
        datadir= fullfile(root,rat,'data');
    end

    if isempty(ps) | ...
        length(ps)<d | isempty(ps{d}) | ...
        length(ps{d})<e | isempty(ps{d}{e}) | ...
        length(ps{d}{e})<t | isempty(ps{d}{e}{t}) 
        error(sprintf('power spectral density was not calculated for tetrode %s [%d %d %d]',rat,d,e,t));
    end

    tstr= sprintf('%s [%d %d %d], itet= %d\n',rat,d,e,t,itet);
    f= ps{d}{e}{t}.freq;
    px= ps{d}{e}{t}.psd;
    dT= diff(ps{d}{e}{t}.win,[],2);

    if cutoff>0
        for itraj=0:3
            tm= getTimes('occ', rat, d, e, itraj, cutoff);
            tmax(itraj+1)= tm(cutoff);
            if ~isfinite(tm(cutoff)); tmax(itraj+1)= max(tm); end
        end
        disp(tmax)

        loadVar(datadir, 'lindistpos', d, rat, 1);
        dt= mean(diff(lindistpos{d}{e}.data(:,1)));
        tind= round((ps{d}{e}{t}.win(:,1)-lindistpos{d}{e}.data(1,1))/dt)+1;
        wintraj= lindistpos{d}{e}.estinfo(tind)+1;
        dT(ps{d}{e}{t}.win(:,1) > tmax(wintraj)')= 0;
%        dT(ps{d}{e}{t}.win(:,1) < tmax(wintraj)')= 0;
    end

    dT= dT/sum(dT);
    ind= find(f>1 & f<30);

    if mtm
        % weighted average of psd's
        pxm= px*dT;
%        pxm= mean(px,2);
        b=fir1(15,.005);
        pxs= filtfilt(b,1,pxm);
%        pxs= pxs/mean(pxs(ind)); %@@
    else
        % filter power spectrum
        b=fir1(200,.01);
        pxs= filtfilt(b,1,px);
    end
    if showPS
        clf;  hold on
        if mtm
%            plot(f(ind), [pxs(ind) pxm(ind)]);
            plot(f(ind), PSfac*[pxs(ind)], 'r');
        else
            plot(f(ind), px(ind));
            plot(f(ind), pxs(ind), 'g');
        end
        xlabel('freq (Hz)');
        ylabel('normalized power');
        axis tight
%        title(tstr);
        if printPS; 
            myprint([1.5 1], sprintf('%s/work/%s-%s-%.2d-%d-%.2d%s',root,figstr,rat,d,e,t,novelStr)); 
        end
        if pausePlot
            disp(tstr);
            pause
        end
    end
    % analyze power spectrum : peak height, width, peak location
    [ix, lx, hy, wx]= anaPeaks(pxs(ind),f(ind), .50, 'right');
    ix= find(isfinite(wx));
%    ix= find(wx<6);
    ix= ix(6<lx(ix) & lx(ix)<10);  % pick peaks in theta band

%     deal with multiple peaks within theta band
    if length(ix)>1; 
        % drop peaks that are less than half as high as the tallest
        mhy= max(hy(ix));
        ix= ix(hy(ix)/mhy > .5);
        nix= length(ix);
        if nix>1
            % merge peaks that are very close and not well separated
            sep= diff(lx(ix));
            for j=1:length(sep)
                i1=ix(j);
                i2=ix(j+1);
                if all(sep(j)<wx([i1,i2]))
                    ix(j)=nan;
                    wx(i2)= wx(i1+1)+sep(j)/2;
                    lx(i2)= mean(lx([i1,i2]));
                    hy(i2)= mean(hy([i1,i2]));
                end
            end

            ix= ix(isfinite(ix));
            if length(ix)>1; 
                warning('Still more than one peak, what now?');
                keyboard
            end

        end
    end
    if isempty(ix); 
        disp(tstr);
        warning('could not find EEG peak within theta band'); 
        ana.theta(itet)= nan;
        ana.max(itet)= nan;
        ana.width(itet)= nan;
        if(mtm)
            ana.stdmax(itet)= nan;
        end
%            keyboard
    else
        ana.theta(itet)= lx(ix);
        ana.max(itet)= PSfac*hy(ix);
        ana.width(itet)= wx(ix);
%        if wx(ix)>6; keyboard; end %@@
        if(mtm)
            ana.stdmax(itet)= PSfac*std(px(ix,:),0,2);
        end
    end
    
    if cutoff>0;
        ana.tmax(itet,:)= tmax;
    end
%    ana.ninvalid(itet)= ps{d}{e}{t}.ninvalid;
    ana.fracnwin(itet)= size(ps{d}{e}{t}.win,1)/ size(ps{d}{e}{t}.winall,1);
    ana.nwin(itet)= size(ps{d}{e}{t}.win,1);
    ana.nwinall(itet)= size(ps{d}{e}{t}.winall,1);
%    keyboard
    ana.Twin(itet)= sum(diff(ps{d}{e}{t}.win,[],2));
    ana.Twinall(itet)= sum(diff(ps{d}{e}{t}.winall,[],2));
    ana.fracTwin(itet)= ana.Twin(itet)/ana.Twinall(itet);

%    keyboard
end

if nargin<4
    if(mtm)
%        anafile= sprintf('%s/work/psana-all-%s%s',root,selectid,novelStr);
%        anafile= sprintf('%s/work/psana-60s-%s%s',root,selectid,novelStr);
        anafile= sprintf('%s/work/psana-x180s-%s%s',root,selectid,novelStr);
    else
        anafile= sprintf('%s/work/psana-pwelch-%s%s',root,selectid,novelStr);
    end
    if ~printPS
        save(anafile, 'ana');
    end
end