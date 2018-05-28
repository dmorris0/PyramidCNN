%segEdgePR
%
%  [posScores, negScores] = segEdgePR( pixEstVals, segLabels, radius )
%
%  Calculate positive and negative scores for boundary estimates from
%  pixEstVals.  
%  pixEstVals: [M x N] estimates for boundary
%  segLabels:  [M x N] labels of segments
%  radius:     [1x1] radius around each true seg boundary to look for best
%                    pixEstVals
%  posScores: [1 x Np] scores for boundary pixels
%  negScores: [1 x Nn] scores for seg interior pixels (excluding
%                      segLabels==0 as well as pixels within radius of
%                      boundary pixels)
%
%  Daniel Morris, Feb 2018
%
%  See also: segPixFeats, getImPatch, segs2edgeInside,
%            calcConPixPatchFilt
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function [posScores, negScores] = segEdgePR( pixEstVals, segLabels, radius )

if numel(pixEstVals) ~= numel(segLabels)
    error('pixEstVals must have same size as segLabels');
end

%first define boundary pixels which we want to estimate:
bpix = segPixFeats( segLabels, 'e', 1, true );

%we want highest pixEstVals for each boundary pixel
%in radius around it
cppFilt = calcConPixPatchFilt( radius, true, false );  
nbpix = applyConPixPatchFilt( single(pixEstVals), cppFilt );
maxVal = max( nbpix, [], 3 );
posScores = maxVal( bpix(:)'==1 );

%find which are negatives pixels
%exclude those within radius of bpix
bkeep = bwdist( bpix==1 ) > radius;
negpix = segLabels > 0 & bkeep;
negScores = pixEstVals( negpix(:)' );


end
