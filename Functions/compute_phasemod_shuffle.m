function out = compute_phasemod_shuffle(numX, phasemod, numSh)

r = rand(numX, numSh);
meanMRVmagSh = nan(numSh,1);
vecangSh = nan(numSh,1);
phasemodSh = nan(numSh,1);
for s = 1:numSh
    rphase = 2*pi*r(:,s);
    meanvec = nanmean(exp(1i*rphase));
    meanMRVmagSh(s,1) = abs(meanvec);
    vecangSh(s,1) = angle(meanvec);
    [~, z] = circ_rtest(rphase);
    phasemodSh(s,1) = log(z);
end
% real-phasemod vs shuff-phasemod distribution
modPctRank = 100*(1-(sum(phasemodSh > phasemod)./numSh));
fprintf('mod > %.02f pct Sh \n', modPctRank)
% output
out.meanMRVmagSh = meanMRVmagSh;
out.vecangSh = vecangSh;
out.phasemodSh = phasemodSh;
out.modPctRank = modPctRank;

end