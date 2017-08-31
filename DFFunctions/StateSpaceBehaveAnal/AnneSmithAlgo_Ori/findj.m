
function [jfinal] = findj(p, n, pcrit)

% Anne Smith Feb 28th 2003

% finds the minumum number of ones in a row you need to
% say with certainty p<pcrit that learning has occurred
% formaula is outlined in Appendix C of Smith et al., 2003

% p       is Probability(correct response)
% pcrit   is the p-value (e.g 0.05)
% n       is the sequence length


res = [];
cnt = 0;

for j = 2:25    % check j-values from 2 to 35

    cnt = cnt +1;

    if(n <= 2*j)   % case of a very short sequence

	f(1) = p^j;
	for i = 2:n-j+1
	    f(i) = p^j*(1-p);
	end
	res = [res; j sum(f)];

    else           % longer sequences

	f(1) = p^j;
	for i = 2:j+1
	    f(i) = p^j*(1-p);
	end
	for i = j+2:n-j+1;
	    xx   = 1:i-j-1;
	    s    = sum(f(xx));
	    f(i) = p^j*(1-p)*(1-s);
	end

	res = [res; j sum(f)];

	clear f
    end 

end

dvalues = find(res(:,2) < pcrit);
jfinal  = dvalues(1) + 1;  % add 1 because j-loop started at 2



