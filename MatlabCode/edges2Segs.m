%edges2Segs
%
%  [labels, pixTree, ilabels] = edges2Segs( density, nonEdgeThresh, minFrac )
%
%  Does segmentation on an edge image.  Builds pixTree for quickshift segmentation
%  density:      [MxNx1] density image with zero for edge and 1 for
%                        non-edge and anything in-between
%  nonEdgeThresh:  non-edge threshold, default 0.95
%                  will have no edges for density > nonEdgeThresh
%  minFrac:        can connect segments to larger segments if their shared
%                  boundary is at least this fraction of the smaller
%                  segment circumference, default 0.1
%  
%  Daniel Morris, Jan 2018
%
%  See also: predProb2Segs, applyConPixPatchFilt, getImPatch, segs2edgeInside,
%            segPixFeats, calcConPixPatchFilt, seg2ConPixPatch,
%            superpixAff, plotEdgeSegs, affinityMergeSuperpix
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function [labels, pixTree, ilabels, fignum] = edges2Segs( density, nonEdgeThresh, minFrac, fignum )

if nargin < 2, nonEdgeThresh = 0.98; end
if nargin < 3, minFrac = 0.1; end
if nargin < 4, fignum = []; end

[ilabels, pixTree] = climbDistim( density, nonEdgeThresh );

labels = mergeSegs( ilabels, density, minFrac, 'larger', fignum );

end

function [pinds, qinds] = connectedPeaks( peaks )

cc = bwconncomp( peaks );
p = cell(1,cc.NumObjects);
q = cell(1,cc.NumObjects);
for its=1:cc.NumObjects
    if numel(cc.PixelIdxList{its}) > 1
        p{its} = cc.PixelIdxList{its}(1)*ones(1,numel(cc.PixelIdxList{its})-1);
        q{its} = cc.PixelIdxList{its}(2:end)';
    end
end
pinds = horzcat(p{:});
qinds = horzcat(q{:});
end

function [labels, pixTree] = climbDistim( density, dthresh )

pixDist = 1.3;

%To avoid connections outside the image, make sure all densities are > 0
offset = -min(density(:))+1;
density = density + offset;
dthresh = dthresh + offset;

%This filter shifts +/- 1 pixel in x and y
[cppFilt,xyoff] = calcConPixPatchFilt( pixDist, false );  %false: only include neighbors, no zero offset

%arrange density of neighbors (along dimension 3):
nbDen = applyConPixPatchFilt( density, cppFilt );

[maxNbDen, maxInd] = max( nbDen, [], 3 );
worseNeighbor = maxNbDen <= density;


%deltaInd: for each shift in x and/or y, how many pixel indices does one shift?
deltaInd = -xyoff(:,2)' - xyoff(:,1)'*size(density,1);
dpind = deltaInd(maxInd);       %find best connection for each pixel
dpind( worseNeighbor ) = 0;    %bad connections mean connect to itsel
[xind,yind] = meshgrid( 1:size(density,2), 1:size(density,1) );
pind = yind + (xind-1) * size(density,1);

pixTree = bsxfun( @plus, pind, dpind );

if dthresh > offset
    %Merge superpixels consisting of adjacent pixels with density > dthresh
    peaks = density > maxNbDen | density >= dthresh;
    [pinds, qinds] = connectedPeaks( peaks );
else
    pinds = [];
    qinds = [];
end

labels = pixTree2Labels( pixTree, pinds, qinds );

end





