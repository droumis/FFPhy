function calcOffset(selectid)
%function calcOffset
%  
%  Determine offset of theta phase angle for calculating correlation coef.

%selectid= 'CA1PE';
xcwin= 0.080; % +/- [msec] window width for xcorr calc

[cl, nc]= collectCellList(selectid);
setRoot;

for ic=1:nc
    rat= cl.rat{ic}; 
    d= cl.cellnum(ic,1); e= cl.cellnum(ic,2); 
    tet= cl.cellnum(ic,3); c= cl.cellnum(ic,4); 

    if ic==1 | ~strcmp(rat, cl.rat{ic-1}) 
        selfile= sprintf('%s/%s/data2/select-%s.mat',root,rat,selectid);
        load(selfile);
    end

    if ic==1 | ~strcmp(rat, cl.rat{ic-1}) | d~= cl.cellnum(ic-1,1)
        fprintf(1, 'calc. phase for %s day %d\n', rat,d);
        bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
        load(bfile);
        sdfile= sprintf('%s/%s/data2/spikedata%.2d.mat',root,rat,d);
        load(sdfile);
    end

    if(0)
        t= behavdata{d}{e}.time;
        dt= mean(diff(behavdata{d}{e}.time));
        nt= length(t);

        mint= t(1);
        maxt= t(end)+dt-1e-10;

        xcn= ceil(xcwin/dt);
        count= zeros(nt, 1);
        spind= find(behavdata{d}{e}.traj(spikedata{d}{e}{tet}{c}.index)>0);
        if isempty(spind); continue; end
        it= floor((spikedata{d}{e}{tet}{c}.time(spind)-mint)/dt);
        it= it(it>0 & it<=nt);
        count(it)= count(it)+1;

        xc= xcorr(count, cos(behavdata{d}{e}.phase+pi)+1, xcn, 'unbiased');

        % so that peak of spiking occurs at pi
        [mxc, ixc]= max(xc);
        dtxc= dt*(ixc-xcn-1);
        %@@
    %    fTheta= ana.theta(ep.EEGrefind(ie));
    %    tTheta= 1/fTheta;
        tTheta= 0.125;
        
        offset= 2*pi*dtxc/tTheta; 
    else
        spind= find(behavdata{d}{e}.traj(spikedata{d}{e}{tet}{c}.index)>0);
        offset= circstat(spikedata{d}{e}{tet}{c}.phase(spind))-pi;
    end


    selid= find(ismember(select.cellnum,[d e tet c],'rows'));
    for ia=1:length(select.a{selid})
        select.a{selid}{ia}.offset= offset;
    end
    selid
    select.a{selid}{:}
    if ic==nc | ~strcmp(rat, cl.rat{ic+1}) 
        save(selfile, 'select');
    end
end
