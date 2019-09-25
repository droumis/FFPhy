function [ac,pac,ir,sr,ns] = idf(z,n)
%function [ac,pac,ir,sr,ns] = idf(z,n)
%   IDF- Identification function
%  IN
%   z = [y u], or [x] (DON'T FORGET TO DTREND Z OR X BEFORE USING)
% where: y - output time series
%        u - input time series
%	 x - univariate series
%        n - number of lags (default=20)
%        
%  OUT
%	 ac  - autocorrelations         -----| UNIVARIATE
%	 pac - partial autocorrelations -----| SERIES
%	 ir  - impulse response
%	 sr  - step response
%	 ns  - noise vector y(t)-ir(t)* u(t)


% -----------------------------------------------------------------
%   By    Edward Katende
%         University of Western Ontario
%	Modified A. Jutan. 2/91 2/94,9/95,1/96
% ------------------------------------------------------------------
	ac=[];pac=[];ir=[];sr=[];ns=[]; % default outputs
 %z=dtrend(z); % DETRENDING DATA RESPONSIBILITY OF USER!
 m=length(z);
 [r,c]=size(z);
 if c>2, error('z cannot have more than 2 columns'), return; end
 
 if nargin<2, n=20; end

				if c==1
	figure(1)				 
      plot(z);
    title('Univariant time series ')

% autocorrelation for n lags starting at lag zero
%
   ac=coxy(z,z,'cor',n); % autocorrelation 
   
	cl= 2*sqrt(1/m); % 95 cl
	figure(2)
		
	subplot(211);
timeax=[0:length(ac)-1];
    stem(timeax,ac),line(timeax,cl*ones(1,length(ac)),'Linestyle','-.')
    line(timeax,-cl*ones(1,length(ac)),'Linestyle','-.')
    line(timeax,-cl*zeros(1,length(ac)),'Linestyle','-')

   title('Autocorrelation ')
   xlabel('lag k  ')

% partial-autocorrelation  for n lags starting at lag one
   ac1(1)=1;
   pac(1)=ac(2);
   for i=2:n-1
     ac1(2:i)=ac(2:i) ;
     p=toeplitz(ac1);
     ac2=ac(2:i+1) ;
     pac1=p\ac2 ;  
     pac(i)=pac1(i);   %partial autocorrelation 
   end;
   pac=[1,pac]'; % add in unity for lag 0 on pac
   	subplot(212);

timeax=[0:length(pac)-1];
    stem(timeax,pac),line(timeax,cl*ones(1,length(pac)),'Linestyle','-.'),
    line(timeax,-cl*ones(1,length(pac)),'Linestyle','-.'),
    line(timeax,-cl*zeros(1,length(pac)),'Linestyle','-')
    
   title('Partial autocorrelation ')
   xlabel('lag k  ')
	subplot(111);
				else
	figure(1)
   if( max(z(:,2))-abs(min(z(:,2))))<0.2; % we have PRBS for input u
     idplot(z)
   else
     subplot(211)
     plot(z(:,1))
     title('Output time series ')
     subplot(212)
     plot(z(:,2))
     title('Input time series ') 
   end

% impulse response
%
	   [ir,R,clir]=cra(z,n,[],0) ; % impulse responses no plot
   
figure(2)
   	clir=clir/2.58*2.0; % 95 CL as opposed to 99 CL
timeax=[0:length(ir)-1];
    stem(timeax,ir),line(timeax,clir*ones(1,length(ir)),'Linestyle','-.')
    line(timeax,-clir*ones(1,length(ir)),'Linestyle','-.')
    line(timeax,-clir*zeros(1,length(ir)),'Linestyle','-')
    
   title(' Impulse response ')
   xlabel('lag k with .95 CL"s ')

% step response
%
figure(3)
   sr=cumsum(ir);
timeax=[0:length(sr)-1];
    stem(timeax,sr),line(timeax,clir*ones(1,length(sr)),'Linestyle','-.'),
    line(timeax,-clir*ones(1,length(sr)),'Linestyle','-.'),
    line(timeax,-clir*zeros(1,length(sr)),'Linestyle','-')
    
   title('Step response ')
   xlabel('lag k with .95 CL"s ')

% ---------------------------------------------------------------------
%    Process noise Identification
% --------------------------------------------------------------------

ns=z(:,1)-filter(ir,1,z(:,2)); % cal noise n(t)= y- ir(t)*u(t)
				end




























