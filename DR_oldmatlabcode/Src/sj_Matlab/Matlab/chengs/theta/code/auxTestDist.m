function auxTestDist(x, sel, opt)
%function auxTestDist(x, sel, opt)
%
% test samples in x{i}
%
%   sel{i}     description of samples
%       label
%
%   opt.test=    {'norm', 'signrank'}
%   opt.plot='cdf'; {'cdf', 'pdf', 'hist'}
%   opt.name        name of variable for internal use
%   opt.title       name of variable w/o special char
%   opt.label       name of variable for labeling in figure
%   opt.ref= 1
%   opt.plotcol= {'k', 'r', ...}
%   opt.size

if nargin<3; opt= []; end
if ~isfield(opt, 'ref'); opt.ref=1; end
if ~isfield(opt, 'plot'); opt.plot='cdf'; end
if ~isfield(opt, 'plotcol');
    plotcol(1:6,:)=[[0 0 0]; [1 0 0]; [0 0 1]; [0 1 0]; [.7 0 1]; [0 1 1]];
else
    plotcol= opt.plotcol;
end
if ~isfield(opt, 'test'); 
    opt.test={'norm', 'signrank', 'ttest'}; 
end

rejstr= {'', ' not'};
nsets= length(x);
%keyboard
figure
maxstr= length(sel{1}.label);
for is=2:nsets; 
    maxstr= max(maxstr, length(sel{is}.label));
end
fprintf(1, '\n%*s  ', maxstr, '');
for it=1:length(opt.test)
    switch opt.test{it}
    case 'norm'
        fprintf(1, 'normal\t');
    case 'signrank'
        fprintf(1, 'signran\t');
    case 'ttest'
        fprintf(1, 'ttest\t');
    otherwise 
        error(['unknown test ' opt.test{it}]);
    end
end
fprintf(1, '\n');

for is=1:nsets
    subplot(ceil(nsets/3), 3, is);
    fprintf(1, '%*s: ', maxstr, sel{is}.label);
    for it=1:length(opt.test)
        switch opt.test{it}
        case 'norm'
            if length(x{is})<4; continue; end
            hp= normplot(x{is});
            set(hp, 'Color', opt.plotcol(is,:));
            [h,p]= lillietest(x{is});
            if isfinite(p); 
                fprintf(1, '%.2g\t', p);
            else
                if h==0; 
                    fprintf(1, '>0.20\t');
                else
                    fprintf(1, '<0.01\t');
                end
            end
        case 'signrank'
            p= signrank(x{is});
            fprintf(1, '%.2g\t', p);
        case 'ttest'
            [h,p]= ttest(x{is});
            fprintf(1, '%.2g\t', p);
        otherwise 
            error(['unknown test ' opt.test{it}]);
        end
    end
    fprintf(1, '\n');
end



