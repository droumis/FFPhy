function [b,ji,sigsq] = arft(y,k,m)
% Usage [b,ji,sigsq] = arft(y,k,m)
% Inputs:
% y = data
% k = maximum order
% m = technique: AIC=1, HQIC=2, BIC=3
% Outputs:
% b = vector of AR coeffs, where AR model is
% Y(t)+b(1)*Y(t-1)+...=e(t)
% ji = estimated order
% sigsq = estimated residual variance
T=length(y);
if m==1
   cr = 2/T;
elseif m==2
      cr=2*log(log(T))/T;
else
   cr=2*log(T)/T;
end
z=fft([y;zeros(size(y))])/T;
z=z.*conj(z);
z=fft(z);
c=real(z(1:(T-1)))/2;
c=c(1:k+1);
b=zeros(k,k);
ac=zeros(k+1,1);
s=zeros(k+1,1);
s(1)=c(1);
ac(1)= log(s(1));
b(1,1)=-c(2)/c(1);
s(2)=s(1)*(1-b(1)*b(1));
ac(2) = ac(1) + log(1-b(1,1)*b(1,1)) + cr;
for i = 1:(k-1)
d=c(i+2)+b(1:i,i)'*c(i+1:-1:2);
e=-d/s(i+1);
b(1:i,i+1)=b(1:i,i)+e*b(i:-1:1,i);
b(i+1,i+1)=e;
s(i+2)=s(i+1)*(1-e*e);
ac(i+2) = ac(i+1) + log(1-e*e) + cr;
end
[max,ji]=min(ac);
ji=ji-1;
if ji > 0
	b=b(1:ji,ji);
else
b=[];
end
sigsq=s(ji+1);
