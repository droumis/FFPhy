function  a=jointcr(x,y,b2,b,func,fac,alpha)
% function a=jointcr(x,y,b2{,b,'func',fac,alpha})
% plots joint confidence regions for two parameters in b2
% These 2  parameters are NOT contained in b
% a= matrix used for contour plots
% b= p-2 vector containing Least Squares estimates of remaining
%     parameters whose joint CR is not required, p=Total # of params
%     If No b dont enter or set to 1
% b2= a 2 vector containg Least Sq estimates of the 2 param whose
%     CR is to be plotted  b2=[ba bb]
% x=  independant variables	used by 'sumss' to calc response y
% y=  n response values i.e. data points at the x settings.
% fac= factor such that plot range=  b2(i) +/- fac*b2(i) i=1,2. {1}
%     If fac contains two #'s they are assumed to be std deviations of params
%     in that case range is b2(i) +/- 5*std(i)
% alpha = F value at which contour is to be plotted. From Statistical
%     Tables such that it corresponds to say 95% Confidence Limits.
%     Look up F{ p,(n-p)} at alpha =.05 One Tail, n=# data points  {0.05}
% func = name of function to evaluate Sum of Squares at given b2,b
%	 of the form q=func(x,y,b2,b). Use 'quotes': default='sumss'
% {}= optional parameters ,{} contains defaults otherwise

% create col vector
y=y(:);[m,n]=size(x);
if m ~= length(y), error(' x,y must have equal # rows '),end;
%
% calc ss contour value f ( Stat Notes p89)
if nargin == 7
    l=m-2;
    amin=feval(func,x,y,[b2(1),b2(2)],b); % calc ss min
    %f=amin.*(1 + 2./l.*(fdistinv(2,l,alpha)));% replaced by stats toolbox
	f=amin.*(1 + 2./l.*(finv(1-alpha,2,l)));
end 
% set up defaults
if nargin <= 6, f=[];end
if nargin<=5 ,fac=1;end
if nargin <=4 ,func='sumss';end
if nargin ==3 ,b=1;end
%
% calc equally spaced plotting values for b2
if length(fac)==1
r1=fac*b2(1);r2=fac*b2(2); % use fractional range
else
r1=5.*fac(1);r2=5.*fac(2);% use range of 5 sigma if sigma's given,
end
b1s=[linspace((b2(1) - r1),(b2(1) + r1),20)];
b2s=[linspace((b2(2) - r2),(b2(2) + r2),20)];
% initialise
k=0;
for ba=b1s
k=k+1;
l=0;
       for bb=b2s;
       l=l+1;
       a(k,l)=feval(func,x,y,[ba,bb],b);
       end
end

% plot contours at 5 heights.if no specific height given
%if nargin ==7  cs=contour(a,[0 f],b1s,b2s);else cs=contour(a,5,b1s,b2s);end
if  nargin ==7 
   cs=contour(b1s,b2s,a,[f f]);
else 
   [cs,h]=contour(b1s,b2s,a,5);
   clabel(cs,h);    % lable all contours

end

title(' CONTOUR PLOT param 1 vs param 2')
xlabel(' param 1');ylabel(' param 2 ');
% add Least Squares point to plot
hold on
plot (b2(1),b2(2),'*')
grid
hold off

