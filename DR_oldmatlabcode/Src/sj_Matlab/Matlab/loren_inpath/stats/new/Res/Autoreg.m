% AUTOREG:  Estimation of the order of an autoregression, using one of three 
%           techniques:
%             AIC, which overestimates with probability about 0.3.
%             HQIC (Quinn), which consistently just estimates correctly.
%             BIC, which is consistent, but could underestimate in small samples.
%           The model is:
%             Y(t)+b(1)*Y(t-1)+...=e(t)
%
%     Usage: [x,b,ord,sigsq] = autoreg(y,maxord,technique)
%
%         y = data vector.
%         maxorder = maximum possible order [default = 10].
%         technique = estimation method:
%                       1 = AIC
%                       2 = HQIC [default]
%                       3 = BIC (most often used)
%         -----------------------------------------------
%         x = filtered ("prewhitened") data vector.
%         b = vector of AR coefficients.
%         ord = estimated order.
%         sigsq = estimated residual variance.
%

% RE Strauss, modified from Barry Quinn's function arft() 

% Hannan, E.J. and Quinn, B.G. (1979). The determination of the order of an
%   autoregression, J. Roy. Statist. Soc. B, 41, 190-195.
% Quinn, B.G. (1980). Order determination for a multivariate autoregression,
%   J. Roy. Statist. Soc. B, 42, 182-185.

function [x,b,ord,sigsq] = autoreg(y,maxord,method)
  leny=length(y);
  if (method==1)
    cr = 2/leny;
  elseif (method==2)
    cr=2*log(log(leny))/leny;
  else
    cr=2*log(leny)/leny;
  end;

  z=fft([y;zeros(size(y))])/leny;
  z=z.*conj(z);
  z=fft(z);
  c=real(z(1:(leny-1)))/2;
  c=c(1:maxord+1);
  b=zeros(maxord,maxord);
  ac=zeros(maxord+1,1);
  s=zeros(maxord+1,1);
  s(1)=c(1);
  ac(1)= log(s(1));
  b(1,1)=-c(2)/c(1);
  s(2)=s(1)*(1-b(1)*b(1));
  ac(2) = ac(1) + log(1-b(1,1)*b(1,1)) + cr;
c
b

  for i = 1:(maxord-1)
i
cval = c(i+1:-1:2)
bval = b(1:i,i)
    d=c(i+2)+b(1:i,i)'*c(i+1:-1:2)';
d
    e=-d/s(i+1);
    b(1:i,i+1)=b(1:i,i)+e*b(i:-1:1,i);
    b(i+1,i+1)=e;
    s(i+2)=s(i+1)*(1-e*e);
    ac(i+2) = ac(i+1) + log(1-e*e) + cr;
  end;

  [max,ord]=min(ac);
  ord=ord-1;
  if (ord > 0)
  	b=b(1:ord,ord);
  else
    b=[];
  end;

  sigsq=s(ord+1);
  x=filter(b,[1],y);

  return;