% rip* ripples
% sci*  sequence compression index

showOutliers= 1;

ripVars= {'Zd', 'COINCFRAC'};
ripLabels= {'z-score', 'co-activity prob.'};
sciInd= [2 5; 3 6; 4 7];
epochLabels= {'famArm', 'novelArm'};
plotcol= {'r', 'g', 'b'};

load SCI_data
for ivar= 1:length(ripVars)
    load(['data_' ripVars{ivar}]);
    for day= 1:3
%        for epoch=1:2 % 1: familiar arm, 2: novel arm
        for epoch=2

            % ripples                SCI
            % IP{set}{arm} 
            ipr= IP{day}{epoch};      ips= AUX{sciInd(day,epoch)}.ip;
            eval(sprintf('ripVar= %s{%d}{%d};', ripVars{ivar}, day, epoch));


            x= dPF{sciInd(day,epoch)}; X= [x, ones(length(x),1)];
            y= lp{sciInd(day,epoch)}';
            xmin= min(x); xmax= max(x);

            b= regress(y, X); 
%            b= [2.1737; -2.7602];
            sresid= y- X*b;
            if showOutliers
                figure(1)
                plot(x,y,'ro', 'markerfacecolor', 'r');
                hold on
                plot([xmin xmax],[xmin 1; xmax 1]*b,'k-');
            end

            ix=0; r= []; s=[]; ip= [];

            for is= 1:length(ips)
                ir=  find(ipr==ips(is));
                if isempty(ir); continue; end
            %    r(is)= COINCFRAC{1}{2}(ir);
                ix= ix+1;
                r(ix)= ripVar(ir);
                s(ix)= sresid(is);
                ip(ix)= ips(is);
            end

            opt.xstr= ripLabels{ivar};
            opt.ystr= 'SCI residual (ms)'; 
            opt.yrange=[-5 130];
            if ivar==1; opt.xrange= [-3 11]; 
            else opt.xrange= [-0.05 0.75]; end
            opt.outname= ['scatter_' ripVars{ivar} '_' epochLabels{epoch} '_day' num2str(day)];
            opt.plot= 'corr'; opt.plotcol= plotcol{day};
            %linrel(r, s);

            if ~showOutliers
                figure
                linrel(r, abs(s), opt);
            end

            if showOutliers
                % "outliers"
                b= [2.1737; -2.7602];
                res= y- X*b;
                figure(2), clf
                hist(abs(res),15)
                [h,b]= hist(abs(res),15)
                %saveas(gcf, 'hist_SCI_resid_outlier');

                rej=find(abs(res)>80);
                valid=find(abs(res)<80);
                fprintf(1, '# of outliers= %d\n', length(rej));
                fprintf(1, 'excluding outliers...\n');

                figure(1); hold on
                plot(x(rej),y(rej),'ko', 'markerfacecolor', 'k');

                xvalid= x(valid); Xvalid= [xvalid, ones(length(xvalid),1)];
                yvalid= y(valid);
                bvalid= regress(yvalid, Xvalid)
                figure(1)
                plot([xmin xmax],[xmin 1; xmax 1]*bvalid,'r-');
                %saveas(gcf, 'seq_comp_pf9_novel-1_outlier');
                keyboard
            end

        end
    end
end


