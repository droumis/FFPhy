% SAMPDIST: Demonstrate sampling distribution of the mean
%
%     Syntax: sampdist(N,iter,distrib)
% 
%           N =       size of random sample drawn from distribution
%           iter =    number of iterations
%           distrib = distribution sampled:
%                       n = normal distribution
%                       u = uniform distribution
%                       b = binary distribution (0,1)
%

function sampdist(N,iter,distrib)
  if (distrib == 'u')                   % Uniform distribution
    xmean = zeros(iter,1);
    for it = 1:iter
      s = unifrnd(0,1,N,1);
      xmean(it) = mean(s);
    end;
    figure(1);
    histgram(xmean);
    putxlab('Mean of X');
    putylab('Relative frequency');

    x = linspace(0,1);
    y = unifpdf(x);
    figure(2);
    plot(x,y);
    axis([0 1 0 1.2]);
    putxlab('X');
    putylab('Relative frequency');

  elseif (distrib == 'n')               % Normal distribution
    xmean = zeros(iter,1);
    for it = 1:iter
      s = randn(N,1);
      xmean(it) = mean(s);
    end;
    figure(1);
    histgram(xmean);
    putxlab('Mean of X');
    putylab('Relative frequency');

    x = linspace(-3.5,3.5);
    y = normpdf(x);
    figure(2);
    plot(x,y);
    axis([-3.5 3.5 0 0.45]);
    putxlab('X');
    putylab('Relative frequency');

  elseif (distrib == 'b')
    xmean = zeros(iter,1);
    for it = 1:iter
      s = round(rand(N,1));
      xmean(it) = mean(s);
    end;
    figure(1);
    histgram(xmean);
    putxlab('Mean of a random sample of X');
    putylab('Relative frequency');

    y = [zeros(10,1);ones(10,1)];
    figure(2);
    histgram(y,[],[],[],'rel');
    putxlab('X');
    putylab('Relative frequency');

  else
    error('Invalid distribution [n,u,b]');
  end;

  return;
