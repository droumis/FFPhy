
%Bayesian implementation of Dynamic Analysis of Learning algorithm in 
%Smith et al., 2004
%run this script within Matlab

%there are possible implmentations 
% flag = 0 if the initial condition is unknown 
% flag = 1 if the initial condition has been estimated previously flag = 1

%data is in a long binary vector I

I = [0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0]'%[0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 1 0 ones(1,20)];
T1 =length(I);
flag = 0

if (flag ==0)

 startp = 0.3;  %initial estimate for starting probability
 T0 = 0;
 filename = 'Modelpart1.txt';
 dataStruct = struct('n', I, 'T1', length(I), 'startp', startp);
 
else
 load xprior
 load xtau
 load Tprior
 startx = xprior;
 T0 = Tprior; %just sets the trial x-axis for plotting purposes
 filename = 'Modelpart2.txt';
 dataStruct = struct('n', I, 'T1', length(I), 'startx', startx, 'tauprior', xtau);
end


%initial guesses for the MC chains%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 chains
clear initStructs
randn('seed',0)
init1 = struct( 'tau', 0.5, 'x', randn(1,length(I))); 
init2 = struct( 'tau', 1,   'x', randn(1,length(I))); 
init3 = struct( 'tau', 15,  'x', randn(1,length(I)));

initStructs(1) =  init1; initStructs(2) =  init2; initStructs(3) =  init3;


%call Winbugs from in matlab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[samples, stats, structArray] = matbugs(dataStruct, ...
		fullfile(pwd, filename), ...
		'init', initStructs, ...
                'nChains', 3, ...
		'view', 1, 'nburnin', 3000, 'nsamples', 1000, ...
		'thin', 10, ...
		'monitorParams', {'p','x','tau','startprior1'}, ...
                'Bugdir', 'C:\Program Files\WinBUGS14');
            
%plotting stuff
subplot(2,1,1)

pdata =[];
for t = 1:length(I)
        sort_samples = sort([samples.p(1,:,t) samples.p(2,:,t) samples.p(3,:,t)]);
        total        = length(sort_samples);
        ll           = sort_samples(fix(0.05*total));  %lower 95%interval
        ml           = sort_samples(fix(0.5*total));
        ul           = sort_samples(fix(0.95*total));
        pdata = [pdata; t ll ml ul];
end


plot(pdata(:,1) + T0, pdata(:,2),'r-'); hold on;
plot(pdata(:,1) + T0, pdata(:,3),'r'); hold on;
plot(pdata(:,1) + T0, pdata(:,4),'r-'); hold on;


%save the values of the latent variable to input to the POST period
xprior     = stats.mean.x(end);
xtau1      = stats.std.x(end);
xtau       = 1/xtau1^2;  %precision of x = 1/variance

if(flag == 1)
allsamplesprior2 = [samples.startprior1(1,:)  samples.startprior1(2,:) samples.startprior1(3,:)];
sortthem = sort(allsamplesprior2);
total    = length(sortthem);
ll       = sortthem(fix(0.05*total));  %lower 95%interval
ml       = sortthem(fix(0.5*total));
ul       = sortthem(fix(0.95*total));

errorbar(T0+0.5,ml,ml-ll,ul-ml); hold on;
plot(T0+0.5, ml ,'+');hold on;
end

save xprior xprior
save xtau xtau
Tprior = T0 + T1;
save Tprior Tprior

%crude MC convergence check
TOOBIGx = find(stats.Rhat.x >1.2);
if(~isempty(TOOBIGx))
    fprintf(2,'not converged \n')
end

    



