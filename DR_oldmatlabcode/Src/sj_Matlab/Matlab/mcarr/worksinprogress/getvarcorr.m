function [out index] = getvarcorr(f, binsize, trjpairs, varargin)

difftetrode = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'difftetrode'
            difftetrode = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%prep variables
for g = 1:length(f(1).epochs) %for each group of epochs
    tmpout{g} = []; index{g} = [];
end

for an = 1:length(f) %for each animal
    for g = 1:length(f(an).epochs) %for each group of epochs
        if (size(f(an).output,2)>=g)
            for e = 1:size(f(an).output{g},2)%for each epoch
                if size(f(an).output{g}(e).lf,2) >= 2
                    
                    %for each cell compute overlap within cell, compute peak
                    for c = 1:size(f(an).output{g}(e).lf,2)
                        trajdata = f(an).output{g}(e).lf{c}.trajdata;
                        
                        %compute peak for all traj
                        alltraj = [];
                        for t = 1:length(trajdata)
                            alltraj = [alltraj; trajdata{t}];
                        end
                        peakrate(c) = max(alltraj(:,5)); %find max of occ normd firing rate
                        totalrate(c) = nansum(alltraj(:,5));
                        
                        %compute overlap between all same turn direction pairs
                        % trjpairs = [1 4; 2 3; 3 5; 4 6; 1 6; 2 5];
                        allovlp = []; allnovlp = [];
                        for t = 1:size(trjpairs,1)
                            trj1 = trajdata{trjpairs(t,1)};
                            trj2 = trajdata{trjpairs(t,2)};
                            if ~isempty(trj1) && ~isempty(trj2)
                                trj1 = trj1(:,5);
                                trj2 = trj2(:,5);
                                ovlp = calcoverlap2(trj1, trj2);   % CORRECT OVERLAP??
                                novlp =  calcoverlap2(trj1, trj2, 'normalize', 1);
                                if ovlp>= 0 && novlp >= 0
                                    allovlp = [allovlp; ovlp];
                                    allnovlp = [allnovlp; novlp];
                                end
                            end
                        end
                        if ~isempty(allnovlp)
                            mnnovlp(c) = nanmean(allnovlp);
                            maxnovlp(c) = max(allnovlp);
                            
                        else
                            mnnovlp(c) = NaN;
                            maxnovlp(c) = NaN;
                        end
                    end
                    
                    %for each pair of cells compute noise corr and overlap
                    %btween cells
                    pairsind = nchoosek(1:size(f(an).output{g}(e).lf,2), 2);
                    
                    % if cells must be on diff tetrode
                    if difftetrode == 1
                        ind = zeros(size(f(an).output{g}(e).lf,2), 4);
                        for c = 1: size(f(an).output{g}(e).lf,2)
                            ind(c,:) = f(an).output{g}(e).lf{c}.index;
                        end
                        difftet = logical(ind(pairsind(:,1),3) - ind(pairsind(:,2),3)); %identify cells pairs with cells on different tetrode
                        pairsind = pairsind(difftet,:); %only include cells pairs with cells on different tetrode
                    end
                    
                    for p =1:size(pairsind,1)%for each pair of cells
                        %cellid
                        c1 = f(an).data{g}{e}(pairsind(p,1),:);
                        c2 = f(an).data{g}{e}(pairsind(p,2),:);
                        %compute noise correlation
                        resid1 = f(an).output{g}(e).resid(pairsind(p,1),:);
                        resid2 = f(an).output{g}(e).resid(pairsind(p,2),:);
                        if sum(~isnan(resid1) & ~isnan(resid2))*binsize>10  % if at least 5sec of values to correlate
                            [rho pval ] = corr(resid1', resid2', 'type', 'Spearman', 'rows', 'complete'); % spearman: no assumptions that residuals are normally distributed
                        else
                            rho = NaN;
                            pval = NaN;
                            
                        end
                        ntimebins = sum(~isnan(resid1) & ~isnan(resid2));
                        
                        %compute overlap btween cells
                        trajdata1 = f(an).output{g}(e).lf{pairsind(p,1)}.trajdata;
                        trajdata2 = f(an).output{g}(e).lf{pairsind(p,2)}.trajdata;
                        alltrjovlp = calcoverlap(trajdata1, trajdata2);
                        alltrjnovlp =  calcoverlap(trajdata1, trajdata2, 'normalize', 1);
                        
                        %compute overlap btween cells per traj
                        allnovlp = []; allovlp = [];
                        for t = 1:max(max(trjpairs)) % for each traj 1-6
                            trj1 = trajdata1{t};
                            trj2 = trajdata2{t};
                            if ~isempty(trj1) && ~isempty(trj2)
                                trj1 = trj1(:,5);
                                trj2 = trj2(:,5);
                                ovlp = calcoverlap2(trj1, trj2);   % CORRECT OVERLAP??
                                novlp =  calcoverlap2(trj1, trj2, 'normalize', 1);
                                if ovlp < 0 && novlp < 0
                                    ovlp = NaN; novlp = NaN;
                                end
                            else
                                ovlp = NaN; novlp = NaN;
                            end
                            allovlp = [allovlp; ovlp];
                            allnovlp = [allnovlp; novlp];
                        end
                        
                        %output
                        cellvalues = [mnnovlp(pairsind(p,1)) mnnovlp(pairsind(p,2)) maxnovlp(pairsind(p,1)) maxnovlp(pairsind(p,2)) peakrate(pairsind(p,1)) peakrate(pairsind(p,2))];
                        tmpout{g} = [tmpout{g}; an f(an).epochs{g}(e,:) pairsind(p,:) rho pval alltrjovlp alltrjnovlp cellvalues allovlp' allnovlp' ntimebins totalrate(pairsind(p,1)) totalrate(pairsind(p,2))];
                        %output: an d e cell# rho pval ovlptweencells novlptweenc mnNovlpwithincell1 mnNovlpwithinc2 maxNovlpwc1 maxNovlpwc2 pk1 pk2 overlaptweencellsTrj1-6  NoverlaptweencellsTrj1-6 numtimebins totalfiring1 totalfiring2
                        %        1  2 3 4  5   6   7            8          9           10                     11        12           13      14   15         16-21        22-27         30       
                        index{g} = [index{g}; an f(an).epochs{g}(e,:) c1 c2];
                        %        an d e tet1 cell1 tet2 cell2
                        %         1 2 3 4      5     6   7
                        
                    end
                end
            end
        end
    end
end

out = tmpout;

end