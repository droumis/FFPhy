%Script to run the subroutines for binomial EM
%Version 1.2
%
%Anne Smith, April 28th, 2003
%
%updated Anne Smith, August 10th, 2004  - adjusted initial variance for startflag=2 case
% 
%variables to be reset by user:
%        mm (1 by number_trials)  vector of number correct at each trial (see below)
%        ll (1 by number_trials)  total number that could be correct at each trial
%        background_prob          probabilty of correct by chance
%        sige                     sqrt(variance) of random walk

%other variables
%        q, s   (vectors)         hidden process and its variance (forward estimate)
%        qnew, signewsq (vectors) hidden process and its variance (backward estimate)
%        newsigsq                 estimate of random walk variance from EM 
%        p      (vectors)         mode of prob correct estimate from forward filter
%        p05,p95   (vectors)      conf limits of prob correct estimate from forward filter
%        b      (vectors)         mode of prob correct estimate from backward filter
%        b05,b95   (vectors)      conf limits of prob correct estimate from backward filter




clear; pack; close all

exampledata = 0;       % flag for example data choice
startflag   = 2;       % 0-to fix initial condition (more likely to give a result)
                       % 1-to estimate initial condition 
                       % 2-to remove xo from likelihood - this means that the latent
					   %   learning process is not started at 0 and allows for alot
					   %   of bias

if(exampledata == 0)   %example of single sequence
 %mm         = [1 0 0 0 0 0 1 0 0 0 1 1 1 1 0 1 0 0 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
 mm          = [0 0 0 0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 1 0 1 ones(1,20)];
 background_prob = 0.5;   
 ll         = ones(size(mm));
else                   %example of proportions
 mm         = [4 3 3 3 2 3 10 11 12 3 5 6 7 8 9 12 14 15 16 18 22 25 26];
 ll         = 30*ones(size(mm));
 background_prob = 0.1;   
end
 I          = [mm; ll];

%starting guess for sige = sqrt(sigma_eps squared)
 sige        = sqrt(0.005) ;                    
 sigsqguess  = sige^2;

%set the value of mu from the chance of correct
 muone           = log(background_prob/(1-background_prob)) ;

%convergence criterion for sigma_eps_squared
 cvgce_crit = 1e-8;

%----------------------------------------------------------------------------------
%loop through EM algorithm 

qguess         = 0;  %starting point for random walk q
number_steps  =  2000;

for jk=1:number_steps

   %do forward filter then backward filter, then EM

   [p, q, s, qold, sold] = recfilter(I, sige, qguess, sigsqguess, muone);

   [qnew, signewsq, a]   = backest(q, qold, s, sold);

   if (startflag == 1)
    qnew(1) = 0.5*qnew(2);   %updates the initial value of the latent process
    signewsq(1) = sige^2;
   elseif(startflag == 0)
    qnew(1) = 0;             %fixes initial value (no bias at all)
    signewsq(1) = sige^2;
   elseif(startflag == 2)
    qnew(1) = qnew(2);       %xo = x1 means no prior chance probability
    signewsq(1) = signewsq(2);
   end


   [newsigsq(jk)]         = em_bino(I, qnew, signewsq, a, muone, startflag);


   qnew1save(jk) = qnew(1);

   %check for convergence
   if(jk>1)
      a1 = abs(newsigsq(jk) - newsigsq(jk-1));
	  a2 = abs(qnew1save(jk) -qnew1save(jk-1));
     if( a1 < cvgce_crit & a2 < cvgce_crit & startflag >= 1)
      fprintf(2, 'EM estimates of RW variance and start point converged after %d steps   \n',  jk)
      break
     elseif ( a1 < cvgce_crit & startflag == 0)
	  fprintf(2, 'EM estimate of RW variance converged after %d steps   \n',  jk)
	  break
     end
   end
 
   sige   = sqrt(newsigsq(jk));
   qguess = qnew(1);
   sigsqguess = signewsq(1);

end


if(jk == number_steps)
 fprintf(2,'failed to converge after %d steps; convergence criterion was %f \n', jk, cvgce_crit)
end

%-----------------------------------------------------------------------------------
%integrate and do change of variables to get confidence limits

[b05, b95, bmid, bmode, pmatrix] = pdistn(qnew, signewsq, muone, background_prob);

%-------------------------------------------------------------------------------------
%find the last point where the 90 interval crosses chance
%for the backward filter (cback)

 cback = find(b05 < background_prob);

 if(~isempty(cback))
  if(cback(end) < size(I,2) )
   cback = cback(end);
  else
   cback = NaN;
  end
 else
   cback = NaN;
 end

%-------------------------------------------------------------------------------------
%plot up the figures

 t=1:size(p,2)-1; 

 figure(1);  clf;

 subplot(211);  
 plot(t, bmode(2:end),'r-');
 hold on;
 plot(t, b05(2:end),'k', t, b95(2:end), 'k');
 if(exampledata == 0)
  hold on; [y, x] = find(mm > 0);
  h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor','k');
  set(h, 'MarkerEdgeColor', 'k');
  hold on; [y, x] = find(mm == 0);
  h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor', [0.75 0.75 0.75]);
  set(h, 'MarkerEdgeColor', 'k');
  axis([1 t(end)  0 1.05]);
 else
  hold on; plot(t, mm./ll,'ko');
  axis([1 t(end)  0 1]);
 end
 line([1 t(end)], [background_prob  background_prob ]);
 title(['IO(0.95) Learning Trial = ' num2str(cback) ' RW variance = ' num2str(sige^2) ]);
 xlabel('Trial Number')
 ylabel('Probability of a Correct Response')



 subplot(212)
 plot(t,1 - pmatrix(2:end),'k')
 line([ 1 t(end)],[0.90 0.90]);
 line([ 1 t(end)],[0.99 0.99]);
 line([ 1 t(end)],[0.95 0.95]);
 axis([1 t(end)  0 1]);
 grid on;
 xlabel('Trial Number')
 ylabel('Certainty')

