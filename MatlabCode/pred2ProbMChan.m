%pred2ProbMChan
%
%  prob = pred2ProbMChan( pred )
%
%  Convert multi-channel binary predictions to probabilities.  Basically
%  soft-max on binary predictions, see mchansoftmax in vl_nnloss
%  pred: [M x N x 2*P x Q]
%  prob: [M x N x P x Q]
%
%  Daniel Morris, Nov 2017
%
%  See also: getTestResults, showNetWeights, showNetPRcurve, runNet,
%            runNetFolder
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function prob = pred2ProbMChan( pred )

n = size(pred,3);
p1 = pred(:,:,1:n/2,:);
p2 = pred(:,:,n/2+(1:n/2),:);
ediff = exp(p1-p2);
prob = ediff ./ (ediff + 1);
%pmax = max(p1,p2);
%prob = exp(p1-pmax) ./ ( exp(p1-pmax) + exp(p2-pmax) );

end

