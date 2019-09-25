function n= getSelectId(num, select)
if nargin== 1
    global fmaux select
    loadFile(fmaux.select);
end
n=1;
ncells= size(select.cellnum,1);
while(n <= ncells & sum(select.cellnum(n,:)==num)~= 4) 
    n=n+1; 
end
if(n > ncells) 
%    num
%    error('could not find requested cell'); 
    n= [];
end

