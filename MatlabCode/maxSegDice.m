%maxSegDice
%
%  [diceVals, segAssoc, precision, recall] = maxSegDice( gtSegs, estSegs )
%
%  Finds association between gtSegs and estSegs that maximizes the total
%  dice score over all gt segs.  Note: this uses Matlab's linear
%  programming function in the optimization toolbox.
%  gtSegs: [M x N] segment labels, 1:Ns segments (zero values imply background are ignored)
%  estSegs: [M x N] estimated segment labels, values between 1 and Ne
%  diceVals: [1xNs] Dice values for each gt segment.  
%            Note: Dice = 2*Precision*Recall/(Precision+Recall)
%  segAssoc: [Ns x 2] first column is list of gt segments, and second
%             column is corresponding est segment
%
%
%  Daniel Morris, Jan 2018
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function [diceVals, segAssoc, precision, recall] = maxSegDice( gtSegs, estSegs )

[diceVals, mseglist, precision, recall] = maxDice( gtSegs, estSegs );

%let's eliminate elements with no gt segment (occurs when gtSegs has
%skipped elements in range 1:Ns) 
tlabs = double(unique(gtSegs));
if tlabs(1)==0, tlabs = tlabs(2:end); end
Ns = numel(tlabs);  %Number of ground truth segments
if numel(diceVals) > Ns
    diceVals = diceVals(tlabs); 
    mseglist = mseglist(tlabs); 
    precision = precision(tlabs); 
    recall = recall(tlabs);
end

segAssoc = [tlabs(:) mseglist(:)];


end

%chooses association between GT and EstSegs that maximizes DICE:
%it uses linear programming to ensure injection between them
function [diceVals, mseglist, precision, recall] = maxDice( gtSegs, estSegs )

Ns = double( max(gtSegs(:)) );
Ne = double( max(estSegs(:)) );
Nmin = min(estSegs(:));
if Nmin < 1, error('estSegs must be between 1 and Ne'); end

%Fill in matching table:
%with first row being the background (label 0):
match = sparse(double(gtSegs(:))+1, double(estSegs(:)), ones(numel(gtSegs),1), Ns+1, Ne );

nest = max( sum(match, 1), 1e-5);  %avoid divide by zero
ngt = max( sum(match, 2), 1e-5); %avoid divide by zero

%No need to match background (gtSegs == 0)
precision = bsxfun(@rdivide, match(2:end,:), nest );  %tp / (tp + fp)
recall = bsxfun(@rdivide, match(2:end,:), ngt(2:end) );      %tp / (tp + fn)

dice = 2 * precision .* recall ./ max( precision + recall, 1e-5 );

%want one seg per GT and each seg should be different:
fprintf('Finding segment association that maximizes dice... ');
[diceVals, mseglist] = maxInjection( dice );  
fprintf(' done\n');

inds = sub2ind([Ns Ne], 1:Ns, mseglist );

diceVals = full(diceVals)';
precision = full(precision(inds))';
recall = full(recall(inds))';

end



