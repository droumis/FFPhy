function out = stack(inputA, inputB)
%out = stack(inputA, inputB)
%Creates a matrix with the two matrices stacked.  If the 2nd
%dimension of the two matrices are equal, this is equivalent to
%[inputA;inputB].  If not, then the smaller matrix is filled in with NaN's
%until the sizes match.

if ~(isnumeric(inputA) & isnumeric(inputB))
   error('Both inputs must be numeric');
end

sizA = size(inputA);
sizB = size(inputB);

if sum(sizA == 0) > 0
  out = inputB;
  return;
elseif sum(sizB == 0) > 0
  out = inputA;
  return;
end

if (length(sizA) ~= length(sizB)) | ...
  sum(sizA == sizB) < (length(sizA) - 1)
  error('Multidimensional stacking confusion.');
end

catInd = find(sizA ~= sizB);
if isempty(catInd) % default to 1
  catInd = 1;
end
out = cat(catInd,inputA,inputB);

return;


try
   out = [inputA;inputB];
catch

   if ((size(inputA,3) > 1) | (size(inputB,3) > 1)) % stack in 3rd dimension!
      if (sizA ~= sizB)
         error('Multidimensional stacking confusion.');
      else
         sizA = size(inputA,1);
         sizB = size(inputB,1);
         if (sizA == sizB)
            out = cat(3,inputA,inputB);
         elseif (sizA > sizB)
            inputB(sizB+1:sizA,:,:) = nan;
            out = cat(3,inputA,inputB);
         else
            inputA(sizA+1:sizB,:,:) = nan;
            out = cat(3,inputA,inputB);
         end
      end
   else
      if (sizA > sizB)
         inputB(:,sizB+1:sizA) = nan;
         out = [inputA;inputB];
      elseif (sizB > sizA)
         inputA(:,sizA+1:sizB) = nan;
         out = [inputA;inputB];
      else
         error('Can not stack');
      end
   end
end

