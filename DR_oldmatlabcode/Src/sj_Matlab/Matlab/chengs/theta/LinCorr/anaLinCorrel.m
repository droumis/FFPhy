%%%key anaLinCorr

% parameters
calcLC=     1;
collectLC=  1;
loadLC =    0;
plotLC=     1;
printout=   0;

if (loadLC)
    load rAll;
else
    prefix={'kyl', 'ter', 'fel', 'sta'};
%    prefix={'kyl', 'ter'};
%    prefix={'kyl'};
    passes=1:16;
    minspikes= 4;

    global fmaux % file manager auxillary variable

    % set search paths
    addpath('/bach/theta/common')
    addpath('/bach/theta/code')
    addpath('/bach/AdaptFilter/main')

    nPasses= length(passes);
    nRats= length(prefix);

    olddir= pwd;
    rAll= cell(nPasses,1);
    for r=1:nRats
        cd(['/bach4/' prefix{r} '/adapt']);
        setLocalOptions;

%    fmaux.selectid= 'all-passtimes-placefields';
    fmaux.selectid= 'all-passtimes-placefields';
    fmaux.select=['select-' fmaux.selectid];

        if (calcLC)
            disp('**** calculating linear correlation coeff ****')
            calcLinCorr
        end
    end
    cd(olddir)

    if collectLC
        disp('**** collecting  linear correlation coeff ****')
        [mAll mIndiv sAll sIndiv]= collectStats(prefix, 'r', passes, minspikes);
        [mdAll mdIndiv dAll dIndiv]= collectStats(prefix, 'delta', passes, minspikes);
%        keyboard
%        fprintf(1, 'pass %d: mean= %f, std= %f, n= %d\n', ...
%            p, mAll.mean(p), mAll.std(p), mAll.n(p));
        data= [passes; mAll.mean(passes)];
        fstring= 'pass %2d: all= %.2f';
        for r=1:nRats
            fstring= [fstring ', ' prefix{r} '= %.2f'];
            data= [data; mIndiv{r}.mean(passes)];
        end
        fstring= [fstring, '\n'];
        fprintf(1, fstring, data);
    
        if collectLC & calcLC
            save rAll data mAll mIndiv sAll sIndiv prefix passes minspikes
        end
    end

end

if plotLC
    figure(1); clf
    subplot(3,2,1);
    errorbar(passes, mAll.mean(passes), mAll.std(passes)./mAll.n(passes));
    title('all');
    axis([passes(1) passes(end) -0.65 -0.15]);
    xlabel('pass');
    ylabel('lin corr coeff');
    for r=1:length(prefix)
        subplot(3,2,2+r);
        errorbar(passes, mIndiv{r}.mean(passes), mIndiv{r}.std(passes)./mIndiv{r}.n(passes));
        title([prefix{r}]);
%        axis([passes(1) passes(end) -0.65 -0.15]);
        axis tight;
        xlabel('pass');
        ylabel('lin corr coeff');
    end
    print('-dps', 'summary');

    % plot histogram of lin corr coeff's
    plotmap=[1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16];
    binwidth= .2;
    binloc=[-1:binwidth:1, inf];
    for r=1:length(prefix)
        figure(1+r); clf
        for i=passes
            if isempty(sIndiv{r}{i}); continue; end
            h=histc(sIndiv{r}{i}, binloc);
            subplot(4,4, plotmap(i));
            bar(binloc(1:end-1)-binwidth/2,h(1:end-1));
            hh = findobj(gca,'Type','patch');
            set(hh,'FaceColor','k','EdgeColor','w')
            axis([-1 1 0 max(h)]);
            if i==1
                title([prefix{r} 'pass ' num2str(i)]);
            else
                title(['pass ' num2str(i)]);
            end
        end
        print('-dps', ['histo-' prefix{r}]);
    end
end

if printout & plotLC
    figs= findobj('Type', 'figure');
    for f=1:length(figs)
        figure(figs(f));
%        orient landscape;
        orient tall;
        print -dps -P811 
    end
end


