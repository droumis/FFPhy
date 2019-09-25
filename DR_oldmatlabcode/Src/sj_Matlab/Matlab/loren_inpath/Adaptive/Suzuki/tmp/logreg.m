function logreg(dataset,n1,n2)

if nargin<2
    n1=1;
    n2=length(dataset);
elseif nargin<3
    n2=n1;
end

ofile=input('Please input filename [default: logreg.xls]:','s');
if isempty(ofile)
    ofile='logreg.xls';
 end
 ofile
fid = fopen(ofile,'w');
fprintf(fid,'filename\t cell\t cond\t stuff\n');
for i = n1:n2
    for j = 1:length(dataset{i}.cell)
%       for k = 1:length(dataset{i}.conds)
       for k = 10:10
          condtrial=dataset{i}.condtrial(k);
            if length(find(dataset{i}.cobj.cond_no==dataset{i}.conds(k))) > 3 & ~isempty(condtrial{1})
                
                logreg=getlogistcorrect(dataset,i,dataset{i}.conds(k));
                fixrate    =dataset{i}.cell{j}.fixrate(dataset{i}.condtrial{k});
                stimrate   =dataset{i}.cell{j}.stimrate(dataset{i}.condtrial{k});
                delayrate  =dataset{i}.cell{j}.delayrate(dataset{i}.condtrial{k});
                resprate   =dataset{i}.cell{j}.resprate(dataset{i}.condtrial{k});
                
                fprintf(fid,'%s\t cell%d\t %d\t LogReg:\t',dataset{i}.name,j,dataset{i}.conds(k));
                %printrow(fid,logreg(:,2), '%.6f\t');
                fprintf(fid,'%s\t cell%d\t %d\t Scene-Base:\t',dataset{i}.name,j,dataset{i}.conds(k));
                %printrow(fid,stimrate, '%.6f\t');
                fprintf(fid,'%s\t cell%d\t %d\t Delay-Base:\t',dataset{i}.name,j,dataset{i}.conds(k));
                %printrow(fid,delayrate, '%.6f\t');
                fprintf(fid,'%s\t cell%d\t %d\t Resp-Base:\t',dataset{i}.name,j,dataset{i}.conds(k));
                %printrow(fid,resprate, '%.6f\t');
                
                stimrateb=stimrate-fixrate;
                delayrateb=delayrate-fixrate;
                resprateb=resprate-fixrate;
                
                mini=min(min([stimrateb,delayrateb,resprateb]));
                maxi=max(max([stimrateb,delayrateb,resprateb]));
                rangi=maxi-mini;
                if mini<0
                    miniadd=mini;
                else
                    miniadd=0;
                end
                logregm=(logreg(:,2)*rangi)+miniadd;
                graphvars=[logregm';stimrateb;delayrateb;resprateb]';
                figure;
                plot(graphvars);
		legend('probcorrect', 'stimrate', 'delayrate', 'resprate');
                axis tight;
                xlabel('trial');
                ylabel('firing rate or % Correct (blue)');
                title(strcat(dataset{i}.name,' cell',int2str(j),' cond ',int2str(dataset{i}.conds(k))));
            else
                fprintf(fid,'%s\t cell%d\t %d\t not enough trials\n',dataset{i}.name,j,dataset{i}.conds(k));
            end
        end
    end
end
fclose(fid);
