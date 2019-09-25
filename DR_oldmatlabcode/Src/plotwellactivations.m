figure;
day=6;
epoch=6;
t3=find(data{1,day}{1,epoch}.Events.Wellevents(:,3)==14|data{1,day}{1,epoch}.Events.Wellevents(:,3)==13);

% well activation from raw dio events

triggers=[13 14];
tpos=[];
for i=triggers+1;
    triggertimes=DIO{1,day}{1,epoch}{1,i}.pulsetimes;
    for j=1:size(triggertimes,1);
        tstart=triggertimes(j,1);
        tend=triggertimes(j,2);
        tstartind=max(find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<=tstart+1000));
        tendind=min(find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000>=tend+1000));
        currtpos=[data{1,day}{1,epoch}.Pos.correcteddata(tstartind:tendind,2) data{1,day}{1,epoch}.Pos.correcteddata(tstartind:tendind,3)];
        tpos=[tpos;currtpos];
        j=j+1;
    end
end
        


plot(data{1,day}{1,epoch}.Pos.correcteddata(:,2),data{1,day}{1,epoch}.Pos.correcteddata(:,3),'g');
xlim([0 160]);
ylim([0 120]);
hold on;

plot(tpos(:,1),tpos(:,2),'.r');

% hold on;
% for i=1:size(t3,1);
%     time=Data{1,day}{1,epoch}.Events.Wellevents(t3(i,1),1);
%     points=find(Data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<=(time));
%     points1=max(points);
%     plot(Data{1,day}{1,epoch}.Pos.correcteddata(points1,2),Data{1,day}{1,epoch}.Pos.correcteddata(points1,3),'.r');
%     hold on;
% end
%     
    
