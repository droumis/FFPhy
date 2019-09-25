function [pl, np, nf]= collectPairs(cl, pairsel, show, expandPF)
%function [pl, np]= collectPairs(cellsel, pairsel)
% Find pairs of cells specified by parisel given a list of cells in cellsel.
%   pairsel
%       all:        all concurrently recorded cells
%       overlap:    overlaping placefields
%       non_overlap:    non-overlaping placefields
%       (diff_tet:   concurrently recorded cells on different tetrodes)
%       (same_tet:   concurrently recorded cells on same tetrodes)
%   n:    # of real pairs
%   pl
%       rat{n}
%       cellnum(n,:)     [num1, num2]
%   for overlaping pf's
%       i{n}(f)       
%       traj{n}(f)    [traj1; traj2]  traj1= row vector
%       pf{n}(:,f)      [x1; y2; x2; y2]

if nargin<2; pairsel= 'all'; end
if nargin<3 | isempty(show); show= 0; end
if nargin<4 | isempty(expandPF); expandPF= 0; end

pl.pairsel= pairsel;
%pl.cellsel= cellsel;

%[cl,nc]=collectCellList(cellsel);
nc= length(cl.rat);
np= 0; nf= 0;
ncp= 0;
for ic=1:nc
    rat= cl.rat{ic};
    num= cl.cellnum(ic,:); d=num(1); e=num(2); tet=num(3); c=num(4);
    traj= cl.traj{ic};
    % find simultaneous recorded cell
    ind= find(strcmp(rat, cl.rat') & ...
        (d== cl.cellnum(:,1)) & (e== cl.cellnum(:,2)) );
    ind= ind(ind>ic);
    if isempty(ind); continue; end
    switch(pairsel)
    case 'all'
        n= length(ind);
        pl.cellnum(np+1:np+n, 1:4)= ones(n,1)*num;
        pl.cellnum(np+1:np+n, 5:8)= cl.cellnum(ind,:);
        [pl.rat{np+1:np+n}]= deal(rat);
        pl.day(np+1:np+n)= cl.day(ind);
        pl.newarm(np+1:np+n)= cl.newarm(ind);
        np= np+length(ind);
        nf= nan;
    case 'overlap'
        if expandPF
            % find overlapping placefields
            for jc=ind'
                pair_overlaps= 0;
                for i=1:length(cl.traj{ic})
                    x= cl.pf{ic}(:,i);
                    for j=1:length(cl.traj{jc})
                        if cl.traj{ic}(i)~= cl.traj{jc}(j); continue; end
                        y= cl.pf{jc}(:,j);
                        if ( (x(1)<y(1) & x(2)<y(1)) | (x(1)>y(2) & x(2)>y(2)) )
                            continue;
                        end
                        if ~pair_overlaps; ncp= ncp+1; end
                        pair_overlaps= 1;

                        nf=nf+1;
                        pl.cellnum(nf, :)= [num cl.cellnum(jc,:)];
                        pl.rat{nf}= rat;
                        pl.day(nf)= cl.day(ic);
                        pl.newarm(nf)= cl.newarm(ic);
                        pl.traj{nf}= cl.traj{ic}(i);
                        pl.pf{nf}= [cl.pf{ic}(:,i);cl.pf{jc}(:,j)];
                        pl.ip(nf)= ncp;
                    end
                end
            end
            np= nf;
        else
            % find overlapping placefields
            for jc=ind'
                nftmp= 0;
                for i=1:length(cl.traj{ic})
                    x= cl.pf{ic}(:,i);
                    for j=1:length(cl.traj{jc})
                        if cl.traj{ic}(i)~= cl.traj{jc}(j); continue; end
                        y= cl.pf{jc}(:,j);
                        if ( (x(1)<y(1) & x(2)<y(1)) | (x(1)>y(2) & x(2)>y(2)) )
                            continue;
                        end
                        nftmp=nftmp+1;
                        pl.cellnum(np+1, :)= [num cl.cellnum(jc,:)];
                        pl.rat{np+1}= rat;
                        pl.day(np+1)= cl.day(ic);
                        pl.newarm(np+1)= cl.newarm(ic);
                        pl.traj{np+1}(nftmp)= cl.traj{ic}(i);
                        pl.pf{np+1}(:,nftmp)= [cl.pf{ic}(:,i);cl.pf{jc}(:,j)];
                    end
                end
                if nftmp>0; np= np+1; end
                nf= nf+nftmp;
            end
        end
%        keyboard
    case 'non_overlap'
        % find cells that have no overlapping placefields
        for jc=ind'
            nftmp= 0;
            overlap= 0;
            for i=1:length(cl.traj{ic})
                x= cl.pf{ic}(:,i);
                for j=1:length(cl.traj{jc})
                    if cl.traj{ic}(i)~= cl.traj{jc}(j); continue; end
                    y= cl.pf{jc}(:,j);
                    if ( (x(1)<y(1) & x(2)<y(1)) | (x(1)>y(2) & x(2)>y(2)) )
                        continue;
                    end
                    overlap= 1;
                end
            end
            if ~overlap
                np= np+1;
                pl.cellnum(np, :)= [num cl.cellnum(jc,:)];
                pl.rat{np}= rat;
                pl.day(np)= cl.day(ic);
                pl.newarm(np)= cl.newarm(ic);
            end
        end
        nf= nan;
    otherwise
        error(['unknown pairsel= ''' pairsel ''''])
    end

end

if show
    for i=1:np
        fprintf(1, '%s [%d %d %d %d], [%d %d %d %d]\n', ...
            pl.rat{i}, pl.cellnum(i,:));
    end
end
