function tetrodenumcells(animdirect,fileprefix)

% used for Annabelle's animals -- counts cells in cellinfo and installs in .numcells field in corresponding tetrode
%  kk 6.12.13

    load([animdirect, fileprefix, 'tetrodeinfo']);
    load([animdirect, fileprefix,'cellinfo']);
    
    o = cellfetch(cellinfo,'numspikes');
    dummy = o.index;
    dayepochtet = unique(dummy(:,1:3),'rows');
    
    for det=1:size(dayepochtet,1)
        
        numcells = 0;
        row = rowfind(dayepochtet(det,:),dummy(:,1:3));   
        
        while row ~=0
            dummy(row,:) = [0 0 0 0];
            numcells = numcells + 1;
            row = rowfind(dayepochtet(det,:),dummy(:,1:3));   
        end 
        
        tetrodeinfo{dayepochtet(det,1)}{dayepochtet(det,2)}{dayepochtet(det,3)}.numcells = numcells;
        
        
    end
        
    disp('done!')
    
    tetinfo = tetrodeinfo;
    
save([animdirect, fileprefix,'tetinfo'], 'tetinfo')