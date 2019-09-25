function testBars(data,trials)

barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);
barsP_mat(data,trials);

% barsP_mat saves the result in some files in this directory, just as BARS
% does (samp_mu)
load -ascii samp_mu;

% put the original data into x,y
x = data(:,1);
y = int32(data(:,2));
yf = data(:,2);

%compute the confidence intervals
confint = zeros(2, size(samp_mu,2));
for ndx = 1:size(samp_mu,2)
    confint(1,ndx) = prctile(samp_mu(:,ndx),1);
    confint(2,ndx) = prctile(samp_mu(:,ndx),99);
end

%get the mean estimate
fmu = mean(samp_mu,1);


%plot the original data
plot(x,y,'.');
hyf = line(x,yf);
set(hyf,'Color','b');
set(hyf,'LineWidth',1);

%plot the estimate
hfmu = line(x,fmu);
set(hfmu,'Color','g');
set(hfmu,'LineWidth',2);

%plot the confidence intervals
hconfl = line(x, confint(1,:));
set(hconfl, 'Color', 'g');
set(hconfl, 'LineStyle', '--');
hconfu = line(x, confint(2,:));
set(hconfu, 'Color', 'g');
set(hconfu, 'LineStyle', '--');
