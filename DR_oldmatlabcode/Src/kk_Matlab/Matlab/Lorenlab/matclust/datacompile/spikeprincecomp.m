function [scores] = spikeprincecomp(waves, ncomp, amps, pass)


covariance = cov(double(reshape(waves(:,:,pass),40*4,[])'));
waves = reshape(waves, 40*4,[]);
coef = pcacov(covariance);
coef = coef(:,1:ncomp);
scores = (double(waves')*(coef));

