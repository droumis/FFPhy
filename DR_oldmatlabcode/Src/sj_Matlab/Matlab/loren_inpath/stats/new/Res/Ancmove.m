% ANCMOVE: Change ancestor function to move root to final position.
%          Rearrange branch-length vector in same way, and then delete branch 
%          length for root.
%
%     Usage: [anc,brlen] = ancmove(anc,brlen)
%

% RE Strauss, 5/23/01

function [anc,brlen] = ancmove(anc,brlen)
  anc = [anc(:)]';
  brlen = [brlen(:)]';

  lenanc = length(anc);
  oldroot = find(anc==0);

  if (lenanc~=length(brlen))
    error('  ANCMOVE: ancestor-fn and branch-length vectors not of equal length.');
  end;
  if (isempty(oldroot))
    error('  ANCMOVE: no root specified in ancestor function.');
  end;

  if (anc(lenanc)~=0)
    newanc = [anc 0];
    i = find(anc==oldroot)
    newanc(i) = (lenanc+1)*ones(1,2)
    newanc(oldroot) = []
    i = find(newanc>oldroot)
    newanc(i) = newanc(i)-1
    anc = newanc;
    brlen(oldroot) = [];
  end;

  if (length(brlen)==lenanc)
    brlen(lenanc) = [];
  end;

  return;
