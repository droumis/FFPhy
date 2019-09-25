% TESTANS: Generate set of randomized answers for a multiple-choice test.  An approximately equal number
%          of each choice is generated by permutation.  Writes results to the screen and to a text file.

% RE Strauss, 7/29/02

clear all;
close all;

letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

disp(' ');
filename = input('Name of output file: ','s');
nquestions = input('Number of test questions: ');
nchoices = input('Number of choices per question: ');

probs = (1/nchoices)*ones(1,nchoices);
max_repeats = 0.07;                                 % Max proportion of repeated answers

resp = 'n';
max_repeats = ceil(max_repeats * nquestions);

while (resp == 'n')
  counts = prbcount(probs,nquestions,[],[],1);
  counts = counts(randperm(length(counts)));
  
  ans = makegrps(1:nchoices,counts);
  ans = ans(randperm(length(ans)));
  
  r = ans(2:end)-ans(1:(end-1));
  repeats = sum(r==0);

  while (repeats > max_repeats)
    counts = counts(randperm(length(counts)));
    ans = makegrps(1:nchoices,counts);
    ans = ans(randperm(length(ans)));
    r = ans(2:end)-ans(1:(end-1));
    repeats = sum(r==0);
  end;
  
  letterans = letters(ans);
  answers = [tostr(1:nquestions), blanks(nquestions)', letterans']
  tofile(answers,filename);
  
  disp(' ');
  resp = input('Satisfactory answers? (y/n)  ','s');
end;

