%pixTree2Labels
%
%  segLabels = pixTree2Labels( pixTree ) 
%  segLabels = pixTree2Labels( pixTree, pinds, qinds ) 
%
%  Partitions image into segments / labels using connected components on
%  pixel tree.
%  pixTree: [MxNxPxQ] each element is index of parent in MxN image.  If
%           points to self then peak of tree.  
%  if P,Q are 1 then can add additional edge connections between pixels.
%  These are specified by pinds, qinds: [1xNc] identifying pixels that
%  should belong in same segments (in addition to the ones defined by
%  pixTree).
%
%  Daniel Morris, Mar 2017
%
%  See also: calcLabels, calcPixScores, calcPixTree
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function segLabels = pixTree2Labels( pixTree, pinds, qinds ) 

[M,N,P,Q] = size(pixTree);

if nargin < 3, pinds = []; qinds = []; end

if P > 1 || Q > 1
    if ~isempty(pinds) 
        error('pinds not valid for P or Q > 1');
    end
    C = mat2cell( pixTree, M, N, ones(1,P), ones(1,Q) );
    sl = cellfun( @(x) pixTree2Labels( x ), C, 'uniformoutput', false );
    segLabels = reshape( cat(3, sl{:} ), [M N P Q] );
    return;
end

pixinds = 1:M*N;
pinds = [pixinds(~isnan(pixTree)) pinds(:)'];
qinds = [pixTree(~isnan(pixTree))' qinds(:)'];
G = sparse(pinds, qinds, 1, M*N, M*N);
[~,segLabels] = graphconncomp(G+G','Directed',false);  %make symmstric
segLabels = reshape(segLabels, M, N);

end
