function auxEEGref(mtm, selectid, varname)
%function auxEEGref(mtm, selectid, varname)
%
%  Determine EEG reference for each epoch as the tetrode with highest peak in
%  EEG power spectrum.


if nargin<1 | isempty(mtm); mtm=1; end
if nargin<2 | isempty(selectid); selectid= 'CA1PE';  end
if nargin<3 | isempty(varname); varname= 'max';  end

setRoot;

load(sprintf('%s/data/tetrodes-%s', root, selectid));
tet= tetrodes;
ntet= length(tet.rat);

epfile= sprintf('%s/data/epochs-%s', root, selectid);
load(epfile);
ep= epochs;
nep= length(ep.rat);

if(mtm)
    anafile= sprintf('%s/work/psana-pmtm-%s',root,selectid);
else
    anafile= sprintf('%s/work/psana-pwelch-%s',root,selectid);
end
load(anafile);


for ie=1:nep
    rat= ep.rat{ie}; d= ep.num(ie,1); e= ep.num(ie,2); 
    ind= find(strcmp(tet.rat', rat) & tet.num(:,1)==d & tet.num(:,2)==e);

    v= ana.(varname)(ind);
    switch varname
    case 'max'
        [vm vi]= max(v);
    case 'width'
        [vm vi]= min(v);
    otherwise
        error('not implemented')
    end
    ep.EEGref(ie,:)= tet.num(ind(vi),:);
    ep.EEGrefind(ie,:)= ind(vi);
    disp(ep.EEGref(ie,3));
end

epochs= ep;
save(epfile, 'epochs');
