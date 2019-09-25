function auxShowVar(var, opts)
%function auxShowVar(var, opts)
%  var.val is matrix of row vectors, one row for each cell/ subject, each column
%           is independent variable
if nargin < 2; opts= []; end
axisopt= 'tight'; 
if(isfield(opts, 'axisopt')) axisopt= opts.axisopt;  end

v= var.val;
nS= size(v,1);
%fh= figure;
%set(fh, 'Position', [82 161 1229 785]);
if(isfield(opts, 'mean') & opts.mean) 
    plot(mean(v));
    ylabel(var.name);
    xlabel('index');
    axis(axisopt);
else 
    if(isfield(opts, 'oneplot') & opts.oneplot) 
        plot(v');
        ylabel(var.name);
        xlabel('index');
        axis(axisopt);
    else
        nx= ceil(sqrt(nS));
        for i=1:nS; 
            subplot(nx,nx,i); 
            plot(v(i,:));
            axis(axisopt);
        end
    end
end

if(isfield(opts, 'title')) 
    title(opts.title);
    th=get(gca, 'Title');
    set(th, 'Interpreter', 'none')
end


