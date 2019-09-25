
function sss_savespksforautocorr2(d,filename)
%%% GET FROM SESSIONFILE

dn=[d '\Autocorr'];
load(filename);

%%%% SITE1

keep=find(session.site1.neurons{2}=='U');
session.site1.spktimes = session.site1.spktimes(keep);
session.site1.neurons{1} = session.site1.neurons{1}(keep);
session.site1.neurons{2} = session.site1.neurons{2}(keep);

session.site1=rmfield(session.site1,{'alldata','alldatastr','thrs','nclu'0'});

for tr=1:length(session.stim);
    xwhist = session.stim{tr}(:,1);
    session.trigtimes(tr)=xwhist(1);
end

savename= ['Sort_autocorr'];
save (savename, 'spkdata');

