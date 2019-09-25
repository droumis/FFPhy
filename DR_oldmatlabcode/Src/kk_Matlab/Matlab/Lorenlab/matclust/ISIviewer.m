function ISIviewer(clustnum,handles)

global clustattrib;
global clustdata;

plotfig = figure('Color',get(0,'DefaultUicontrolBackgroundColor'), ... 
        'Tag','viewISIfig','NumberTitle','off','Name',['ISI viewer: Cluster ',num2str(clustnum)]); 
    
h1 = axes('XScale','log');
hold on
index = clustattrib.clusters{clustnum}.index;
times = clustdata.params(index,1);
times = times/(clustdata.UnitsPerSec);
ISI = diff(times);
ISI = ISI(find(ISI<.5));
edges = [[.001:.0001:.5]];
N = histc(ISI,edges);


if (sum(N))
    phandle = bar(edges,N,'histc');
    set(phandle,'LineStyle','none');
    xlabel('Time between events (seconds)');
    ylabel('Number of events');
    %set(phandle,'XScale','log');
end