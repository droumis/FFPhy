% runavechange.m --- computes behavioral and neural changepoints
% makes running average of neural data (scene and delay mimus baseline)with a 5 trials window
%
%	cfile - name of cortex object file
%	n     - number of bins to use for learning curve
%	b_dur - baseline duration
%	s_dur - scene duration
%	d_dur - delay duration
%	level - level of significance 
%  	ofile - output file 
% baseline= -300, 0, 29
%
%   batin and norefs is only used by the bat file
%   batin is the file it writes to-sent by the bat file
%   norefs is whether or not to use the reference scenes- 1 is don't
%   0 is do
%works: checked may 2001, sw. 
% output in work. 

function runavechange(cfile,ofile,batin,norefs)

b_dur = 300;
s_dur = 500; 
d_dur = 700;

n=1;
if nargin<3
%dfile = 
%c = cortex(dfile,[1,2]);
%c = cortex('C110600.4',1); 
c = cortex('C110600.4',1); 

norefs=0;
end
if nargin==1
%   ofile=strcat('runchnge',cfile(1:7));
end
    ofile = 'runchange';
cond_list = sort(unique(c.cond_no));
file = strcat(ofile,'.ind');

if nargin>=3
   dfile=cfile;
   fp=batin;
   if nargin==3
      norefs=0
   end
else
	fp = fopen(file,'w');
end
goodcond_no=[];
if norefs
   for i=1:length(cond_list)
      if cond_list(i)>70
         goodcond_no=[goodcond_no;cond_list(i)];
      end
   end
   cond_list=goodcond_no;
end

[row,col] = find(c.codes == 105);

fno=3;
counting=0;

for i = 1:length(cond_list)
      tvec = find((c.cond_no == cond_list(i)) & ((c.response == 0 ) | (c.response == 6)));
   %tvec = subtract(tvec,intersect(row,tvec));
   
   tvec = setdiff(tvec,intersect(row,tvec));


if(~isempty(tvec))	
	r = c.response(tvec);
	amt=find((r==0) | (r==6));
      % list will be an array with 1's for correct and 0's for wrong
      list=r(amt)==0;
   new_list=r(amt);   % r and new_list are the same!
   binwidth=5;
   run_avg=zeros(length(list)-binwidth+1,1);

	for q=1:length(list)-binwidth+1 %for the length of run_avg...
   	run_avg(q)=mean(list(q:q+binwidth-1));  
   end

   

	r_base  = rules('',-b_dur,0,29,cond_list(i),-1,-1,tvec);
	r_scene = rules('',0,s_dur,29,cond_list(i),-1,-1,tvec); 
   r_delay = rules('',0,d_dur,30,cond_list(i),-1,-1,tvec);

	% find discrete behavioral change point and learning curve for each condition

 
		
	len = length(trialfiringrate(c,1,b_dur,r_base));
	
	if((~isempty(find(c.clusters == 1))) & (len > n))
      
		b = trialfiringrate(c,1,b_dur,r_base);
		s = trialfiringrate(c,1,s_dur,r_scene);	
		d = trialfiringrate(c,1,d_dur,r_delay);	
      
		b_1 = avg_compress(trialfiringrate(c,1,b_dur,r_base), n);
		s_1 = avg_compress(trialfiringrate(c,1,s_dur,r_scene),n);	
		d_1 = avg_compress(trialfiringrate(c,1,d_dur,r_delay),n);	
      
      binwidth=5;
   	b_1ra=zeros(length(b)-binwidth+1,1);
		for q=1:length(b)-binwidth+1 %for the length of run_avg...
   		b_1ra(q)=mean(b(q:q+binwidth-1));  
   	end
      
   	s_1ra=zeros(length(s)-binwidth+1,1);
		for q=1:length(s)-binwidth+1 %for the length of run_avg...
   		s_1ra(q)=mean(s(q:q+binwidth-1));  
   	end
      
   	d_1ra=zeros(length(d)-binwidth+1,1);
		for q=1:length(d)-binwidth+1 %for the length of run_avg...
   		d_1ra(q)=mean(d(q:q+binwidth-1));  
   	end
      
      difscb1 = (s_1ra - b_1ra);
      difdelb1 = (d_1ra - b_1ra);
      
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t', 'run_avg');
      printrow(fp,run_avg, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t', 'base1');
      printrow(fp, b_1, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t','base1 run_avg');
      printrow(fp,b_1ra,'%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t', 'scene1');
      printrow(fp, s_1, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t','scene1 run_avg');
      printrow(fp,s_1ra,'%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t', 'delay1');
      printrow(fp, d_1, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t','delay1 run_avg');
      printrow(fp,d_1ra,'%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t', 'difscenebase');
      printrow(fp, difscb1, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 1', cond_list(i));
      fprintf(fp,'%s\t', 'dfdelbase');
      printrow(fp, difdelb1, '%.2f\t');
      
      x1=[ difscb1, difdelb1];
      if nargin<3
      if (cond_list(i) < 71) & (cond_list(i) > 65)
         q = cond_list(i)-66; %we're only dealing with 67-70, so this shouldn't be a problem
         figure(1);
         subplot(2,4,q);
         plot(run_avg);
         axis tight;
      	xlabel('trial');
      	ylabel('% Correct');
      	title(strcat(c.name,'----cell 1----',int2str(cond_list(i))));
         subplot(2,4,q+4);
         plot(x1);
         xlabel('trial');
      	ylabel('firing rate diff');
      	title(strcat(c.name,'----cell 1----',int2str(cond_list(i))));
      end
      if (cond_list(i)>71)
         counting=counting+1;
         if counting > 3
            counting = 1;
            fno=fno+2;  %because of the second channel- add 2 instead of 1 (later fno+1 is used)
         end
			figure(fno);   %generally 3, possibly 5 if there are more than three new scenes
         subplot (2, 3, counting);
         plot (run_avg);
         xlabel('trial');
      	ylabel('% Correct');
      	title(strcat(c.name,'----cell 1----',int2str(cond_list(i))));
         subplot(2,3, counting+3);
         plot(x1);
         xlabel('trial');
      	ylabel('firing rate diff');
      	title(strcat(c.name,'----cell 1----',int2str(cond_list(i))));
      end
   end  %for "if nargin<3"   
      
      
end
	
	len = length(trialfiringrate(c,2,b_dur,r_base));
	
	if(~isempty(find(c.clusters == 2)) & (len > n)  )
		b2 = trialfiringrate(c,2,b_dur,r_base);
		s2 = trialfiringrate(c,2,s_dur,r_scene);	
		d2 = trialfiringrate(c,2,d_dur,r_delay);	
      
		b_2 = avg_compress(trialfiringrate(c,2,b_dur,r_base), n);
		s_2 = avg_compress(trialfiringrate(c,2,s_dur,r_scene),n);	
		d_2 = avg_compress(trialfiringrate(c,2,d_dur,r_delay),n);	
      
      binwidth=5;
   	b_2ra=zeros(length(b2)-binwidth+1,1);
		for q=1:length(b2)-binwidth+1 %for the length of run_avg...
   		b_2ra(q)=mean(b2(q:q+binwidth-1));  
   	end
      
   	s_2ra=zeros(length(s2)-binwidth+1,1);
		for q=1:length(s2)-binwidth+1 %for the length of run_avg...
   		s_2ra(q)=mean(s2(q:q+binwidth-1));  
   	end
      
   	d_2ra=zeros(length(d2)-binwidth+1,1);
		for q=1:length(d2)-binwidth+1 %for the length of run_avg...
   		d_2ra(q)=mean(d2(q:q+binwidth-1));  
   	end
      difscb2 = (s_2ra - b_2ra);
      difdelb2 = (d_2ra - b_2ra);

      
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t', 'run_avg');
      printrow(fp,run_avg, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t', 'base2');
      printrow(fp, b_2, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t','base2 run_avg');
      printrow(fp,b_2ra,'%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t', 'scene2');
      printrow(fp, s_2, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t','scene2 run_avg');
      printrow(fp,s_2ra,'%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t', 'delay2');
      printrow(fp, d_2, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t','delay2 run_avg');
      printrow(fp,d_2ra,'%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t', 'difscenebase2');
      printrow(fp, difscb2, '%.2f\t');
      fprintf(fp,'%s\t %s\t %d\t',getfname(cfile),'chan 2', cond_list(i));
      fprintf(fp,'%s\t', 'dfdelbase2');
      printrow(fp, difdelb2, '%.2f\t');
          
          x1=[ difscb2, difdelb2];
      if nargin<3
      if (cond_list(i) < 71) & (cond_list(i) > 65)
         q = cond_list(i)-66;
         figure(2);
         subplot(2,4,q);
         plot(run_avg);
         axis tight;
      	xlabel('trial');
      	ylabel('% Correct');
      	title(strcat(c.name,'----cell 2----',int2str(cond_list(i))));
         subplot(2,4,q+4);
         plot(x1);
         xlabel('trial');
      	ylabel('firing rate diff');
      	title(strcat(c.name,'----cell 2----',int2str(cond_list(i))));
      end
      if (cond_list(i)>71)
         figure(fno+1);
         subplot (2, 3, counting);
         plot (run_avg);
         xlabel('trial');
      	ylabel('% Correct');
      	title(strcat(c.name,'----cell 2----',int2str(cond_list(i))));
         subplot(2,3, counting+3);
         plot(x1);
         xlabel('trial');
      	ylabel('firing rate diff');
      	title(strcat(c.name,'----cell 2----',int2str(cond_list(i))));
      end
    end %for if nargin<3  
      
	  	end




       
      

  end
end

if nargin<3
fclose(fp);
end
