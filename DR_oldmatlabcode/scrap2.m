
dayfilter = [1:8];





for i = length(animals);
    for k = 1:length(pfs{1}(i).output{1}(1,:));
        currind = [];
        if length(pfs{1}(i).output{1}(1,k).trajdata) ==4 && length(pfs{2}(i).output{1}(1,k).trajdata) ==4;  % if all trajs exist for both mod/unmod
            
            %mod gather trajs
            currind = pfs{1}(i).output{1}(1,k).index;
            
            outleftmod = pfs{1}(i).output{1}(1,k).trajdata{1,1}(:,5);
            outrightmod= pfs{1}(i).output{1}(1,k).trajdata{1,3}(:,5);
            
            inleftmod = pfs{1}(i).output{1}(1,k).trajdata{1,2}(:,5);
            inrightmod = pfs{1}(i).output{1}(1,k).trajdata{1,4}(:,5);
            
            %compute corrcoef
            outmodmintrunc = min(length(outleftmod), length(outrightmod));
            [outmodR outmodP] = corrcoef(outrightmod(1:outmodmintrunc),outleftmod(1:outmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
%             outmodAllR(k,5)  = outmodR(2,1);
%             outmodAllR(k,[1 3 4]) = currind(1,[1 3 4]);
           AllR{1}(k,5)  = outmodR(2,1);
            AllR{1}(k,[1 3 4]) = currind(1,[1 3 4]);
            outmodAllP(k,1) = outmodP(2,1);
            
            inmodmintrunc = min(length(inleftmod), length(inrightmod));
            [inmodR inmodP] = corrcoef(inrightmod(1:inmodmintrunc),inleftmod(1:inmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
%             inmodAllR(k,5)  = inmodR(2,1);
%             inmodAllR(k,[1 3 4]) = currind(1,[1 3 4]);
            AllR{2}(k,5)  = inmodR(2,1);
            AllR{2}(k,[1 3 4]) = currind(1,[1 3 4]);
            inmodAllP(k,1) = inmodP(2,1);
            
            %unmod
            outleftunmod = pfs{2}(i).output{1}(1,k).trajdata{1,1}(:,5);
            outrightunmod= pfs{2}(i).output{1}(1,k).trajdata{1,3}(:,5);
            
            inleftunmod = pfs{2}(i).output{1}(1,k).trajdata{1,2}(:,5);
            inrightunmod = pfs{2}(i).output{1}(1,k).trajdata{1,4}(:,5);
            
            %compute corrcoef
            outunmodmintrunc = min(length(outleftunmod), length(outrightunmod));
            [outunmodR outunmodP] = corrcoef(outrightunmod(1:outunmodmintrunc),outleftunmod(1:outunmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
%             outunmodAllR(k,5)  = outunmodR(2,1);
%             outunmodAllR(k,[1 3 4]) = currind(1,[1 3 4]);
            AllR{3}(k,5)  = outunmodR(2,1);
            AllR{3}(k,[1 3 4]) = currind(1,[1 3 4]);
            outunmodAllP(k,1) = outunmodP(2,1);
            
            inunmodmintrunc = min(length(inleftunmod), length(inrightunmod));
            [inunmodR inunmodP] = corrcoef(inrightunmod(1:inunmodmintrunc),inleftunmod(1:inunmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
%             inunmodAllR(k,5)  = inunmodR(2,1);
%             inunmodAllR(k,[1 3 4]) = currind(1,[1 3 4]);
            AllR{4}(k,5)  = inunmodR(2,1);
            AllR{4}(k,[1 3 4]) = currind(1,[1 3 4]);
            inunmodAllP(k,1) = inunmodP(2,1);
            
            %append index and collapse across epochs within day
            
            
        else
            pfs{1}(i).output{1}(1,k).index %print any cells skipped
            
        end
    end
    
    %day chunking
    
                    for mu = 1:length(AllR);
                    clear  tmpnum x y z c ia ic gu
                    currindexmod=[]; cntallcells = 0; total_fracunder3 = []; pu =0;
                    nidxs_an = length(AllR{mu});
                    tmpnum = nan(length(AllR{mu}(:,1)),1);
                    for hu = 1:length(AllR{mu}(:,1));
                        x = AllR{mu}(hu,[1 3 4]);
                        y = sprintf('%d',x);
                        z = str2num(y);
                        tmpnum(hu) = z;
                    end
                    [c ia ic] = unique(tmpnum, 'stable');
                    for gu = 1:length(ia);
                        corrdata_allep{mu}(gu,5) = nanmean(AllR{mu}((find(ic == gu)),5));
                        corrdata_allep{mu}(gu,[1 3 4]) = AllR{mu}((ia(gu)),[1 3 4]);
                    end
                    end

    
  
    
    figure
    subplot(2,3,1);
    hold on
    bar(1, mean(AllR{1}(:,5)), 'b', 'EdgeColor', 'none')
    bar(2, mean(AllR{3}(:,5)),'r', 'EdgeColor', 'none')
    errorbar2([1 2], [mean(AllR{1}(:,5)) mean(AllR{3}(:,5))],  [stderr(AllR{1}(:,5)) stderr(AllR{3}(:,5))] , 0.3, 'k')
    xlim([0.3 2.7])
    set(gca, 'fontsize', 24)
    set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
    ylabel('corrcoef R')
    [h p] = kstest2(AllR{1}(~isnan(AllR{1}(:,5)),5),AllR{3}(~isnan(AllR{3}(:,5)),5));
    rp = ranksum(AllR{1}(:,5),AllR{3}(:,5));
    title({['CorrCoef Outbound' animals{1,i}]; sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p)})
    
    subplot(2,3,4);
    hold on
    bar(1, mean(AllR{2}(:,5)), 'b', 'EdgeColor', 'none')
    bar(2, mean(AllR{4}(:,5)),'r', 'EdgeColor', 'none')
    errorbar2([1 2], [mean(AllR{2}(:,5)) mean(AllR{4}(:,5))],  [stderr(AllR{2}(:,5)) stderr(AllR{4}(:,5))] , 0.3, 'k')
    xlim([0.3 2.7])
    set(gca, 'fontsize', 24)
    set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
    ylabel('corrcoef R')
    [h p] = kstest2(AllR{2}(~isnan(AllR{2}(:,5)),5),AllR{4}(~isnan(AllR{4}(:,5)),5));
    rp = ranksum(AllR{2}(:,5),AllR{4}(:,5));
    title({['CorrCoef Inbound__' animals{1,i}]; sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p)})
    
    
                    %plot corrcoef over days
                for j = [str2num(dayfilter)];
                    daymod{j} = sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5);
                    dayunmod{j} = sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5);
                    if runnadal == 1;
                        if j == 8; %nadal
                            meanday = [mean(sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5)) mean(sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5))];
                        else
                            meanday = [meanday; [mean(sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5)) mean(sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5))]];
                        end
                    else
                        if j == 1; %HP
                            meanday = [mean(sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5)) mean(sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5))];
                        else
                            meanday = [meanday; [mean(sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5)) mean(sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5))]];
                        end
                    end
                    
                end
                
                subplot(2,2,3);
                hold on
                bar([str2num(dayfilter)], meanday, 'EdgeColor','none') %hpa animals nadal
                %                 errorbar2(1:ndays, meanday,  [stderr(mod) stderr(unmod)] , 0.3, 'k') %add error bars later if something interesting
                if runnadal == 1;
                    xlim([7.5 17.5]) %nadal
                else
                    xlim([0 8.5]) %HP
                end
                %
                set(gca, 'fontsize', 24)
                
                set(gca, 'xtick',[str2num(dayfilter)]) %nadal
                %                 set(gca, 'xtick', 1:ndays)
                ylabel('Proportion field size')
                xlabel('days')
                title('Proportion field size >3Hz over days');
    
    
end




