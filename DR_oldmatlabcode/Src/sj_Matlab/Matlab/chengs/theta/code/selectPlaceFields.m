function selectPlaceFields(outfile)
% function selectPlaceFields 
% 
% Select individual place fields, even if one cell has multiple.
% Here, the "place field" is defined as position interval with the following
% properties:

% 1. The peak average firing rate exceeds minpeakrate
% 2. The average firing rate at every point within the place field exceeds
% (edgeratio) x (peak rate) or maxEdgeRate, whichever is less
% 3. The average velocity at every point withing the place field exceeds minvel 
% 4. The width of the place field is at least minwidth and at most maxwidth.

debugging= 0;

%minpeakrate= 10; % minimum firing rate [Hz]
minpeakrate= 3; % minimum firing rate [Hz]
edgeratio= .10;
maxEdgeRate= 20;
%minvel= 10;  % minimum velocity [cm/s]
minvel= 2;  % minimum average velocity [cm/s]
minwidth= 1; % [cm]
%maxwidth= 50; % [cm]
maxwidth= 170; % [cm]  % entire track

traj=[1, 2; 2 1; 1 3; 3 1];
global fmaux linmeanvel behavdata nbinx select task tetloc 
global linbehav spikes rate
loadVar(fmaux.datadir, 'tetloc', 0, fmaux.prefix, 1);
if debugging
    load(fmaux.select);
    fmaux.nCells=size(select.cellnum,1);
else
    loadFile(fmaux.select);
    resetCellList;
end

if ~isfield(select, 'x');
    select.x= {};
end
n= 1;

nPF= [];
while 1
    [d,e,t,c]= getNextCell;
    if isempty(d) | isempty(e) | isempty(t) | isempty(c); break; end
    n= fmaux.currentCell;

%    loadVar(fmaux.datadir, 'spikes', d, fmaux.prefix, 1);
    
    if n>length(select.a) | isempty(select.a{n}) 
        continue;
    end
    nana= length(select.a{n});
    if n>length(select.x) select.x{n}= cell(1,nana); end

    % read in position, rate, and velocity information
    loadVar(fmaux.data2dir, 'rate', 0);
    loadVar(fmaux.data2dir, 'linmeanvel', 0);
    loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, true);

    pixelsize= task{d}{e}.pixelsize;
    
    % save selection criteria and independent variable, so that they can be
    % duplicated (in case there's more than one place field per traj)
    olda= select.a{n};
    oldx= select.x{n};
    select.a{n}= {};
    select.x{n}= {};
    nPFcell= 0;

    rdata= rate{d}{e}{t}{c}.data;
    linpos= rdata(:,1);

    for ia=1:length(olda) % loop over trajectories
        imindone=[];
        imaxdone=[];
        nPFtraj= 0;
        it= olda{ia}.traj;

        ratetmp= rdata(:,it+2);

       %@@ Loren's old rate calc
    %    loadVar(fmaux.datadir, 'linbehav', d, fmaux.prefix, true);
%       linpos= linbehav{d}{e}{t}{c}.data{traj(it+1,1),traj(it+1,2)}(:,1)*pixelsize;
%       rate= linbehav{d}{e}{t}{c}.data{traj(it+1,1),traj(it+1,2)}(:,2);
                                

        % convert pixels/s to cm/s;  some vel could be nan, set those to 0
        vel= linmeanvel{d}{e}.data{it+1}(:,2)*pixelsize;
        vel(find(isnan(vel)))= 0;

        %%%% find place fields %%%%

        % interpolate rate and velocity with cardinal splines
        x= linpos(1):.1:linpos(end);
        r= interp1(linpos,ratetmp,x,'spline');
        tmpvel= linmeanvel{d}{e}.data{it+1}(:,1)*pixelsize;
        v= interp1(tmpvel,vel,x,'spline');
        nx= length(x);

        if (debugging) % debugging output
            global spikedata behavdata
            loadVar(fmaux.data2dir, 'behavdata', d);
            loadVar(fmaux.data2dir, 'spikedata', d);
            sd= spikedata{d}{e}{t}{c};
            bd= behavdata{d}{e};

            disp('==============');
            figure(1)
            subplot(2,1,2);
            trajs= bd.traj(sd.index);
            tind= find(trajs== it);
            plot(sd.linpos(tind), bd.phase(sd.index(tind)), '.'); 
            axis([0 160, 0 2*pi]);
            subplot(2,1,1);
            tstring= sprintf('[%d %d %d %d], cellnum= %d, traj= %d\n',d,e,t,c,n,it);
            plot(x,[r;v]); 
            axis([0 160, 0 max(max(r),max(v))]);
            legend({'rate', 'vel'});
            title(tstring);
            pause
        end

        % find peaks of firing rate
        ipeaks= findPeaks(r);
        % sort peaks with highest first
        [tmp isort]= sort(r(ipeaks));
        ipeaks=ipeaks(fliplr(isort));

        % to keep track of which parts of the trajectory have already been
        % covered, we set the rate within them to -1
        for ip= ipeaks
            if debugging; fprintf(1,'%.1f: ', x(ip)); end
%            if r(ip) < 0; 
%                fprintf(1,'already covered');
%            end
            % check rate and vel at peak
            if r(ip) < minpeakrate; 
                if debugging; fprintf(1,'minpeakrate\n'); end
                continue; 
            end

            % find edges of putative place fields
            ratecut= edgeratio*r(ip);
            if ratecut > maxEdgeRate;
                ratecut= maxEdgeRate;
            end
            itmp= find(r < ratecut);
            imin= max(itmp(find(itmp < ip)));
            if isempty(imin); 
                imin=1; 
            else
                imin= imin+1;
            end
            imax= min(itmp(find(itmp > ip)));
            if isempty(imax); 
                imax=length(r); 
            else
                imax= imax-1;
            end
            r(imin:imax)= -1;

            minx= x(imin); maxx= x(imax);
            width= maxx-minx;
            % check width
            if width < minwidth; 
                if debugging; fprintf(1,'minwidth\n');end
                continue; 
            end
                
            if find(v(imin:imax) < minvel); 
                if debugging; fprintf(1,'minvel\n');end
                continue; 
            end
            if width > maxwidth;  
                if debugging; fprintf(1,'maxwidth\n');end
                continue; 
            end

            if debugging
                hold on
                plot(minx,0,'sr'); 
                plot(maxx,0,'or');
                fprintf(1,'place field [%.1f %.1f]!\n', minx, maxx);
                hold off
            pause
            end

            % accept place field and save
            nPFtraj= nPFtraj+1;
            nPFcell= nPFcell+1;
            select.x{n}{nPFcell}= oldx{ia};
            select.a{n}{nPFcell}= olda{ia};
            select.a{n}{nPFcell}.linpos= [minx; maxx];
%            select.a{n}{nPFcell}.phase= [0; 2*pi];
        end % ip (peaks)
        nPF(n,ia)= nPFtraj;
    end % ia (= traj)
end

[select,totalana]= cleanSelect(select);
% save to disk
save('select-pf10', 'select','totalana')

