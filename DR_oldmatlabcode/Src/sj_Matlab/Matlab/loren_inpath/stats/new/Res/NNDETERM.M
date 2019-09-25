% NNDETERM: For a fully connected feedforward neural network, determines 
%           either the number of training vectors needed (for a given 
%           determination ratio) or the determination ratio (given the 
%           number of training vectors available).
%
%     Usage: function x = nndeterm(input,hidden,output,trainvect,determ)
%
%           input =     number of input nodes
%           hidden =    vector of numbers of nodes in hidden layers
%           output =    number of output nodes
%           trainvect = number of training vectors
%           determ =    determination ratio (number of training vectors / exact 
%                         number needed for an exactly determined network
%
%           Either trainvect or determ must be provided, and the other (the 
%           value of which will be returned) set equal to zero.

%     Reference: Carpenter,W.C. & M.E. Hoffman, "Training backprop neural;
%                   networks", AI Expert 10(3):30-33.

% RE Strauss, 2/8/95

function x = nndeterm(input,hidden,output,trainvect,determ)
  if (nargin<5)
    error('  NNDETERM: All five arguments required');
  end;

  nlayers = length(hidden);
  nparam = hidden(1) * (input+1);
  if (nlayers > 1)
    for i=2:nlayers
      nparam = nparam + (hidden(i) * (hidden(i-1)));
    end;
  end;
  nparam = nparam + (output * (hidden(nlayers)+1));

  if (trainvect == 0)
    x = determ * nparam;
  else
    x = trainvect / nparam;
  end;
    
  return;
