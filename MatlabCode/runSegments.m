%runSegments
%
%  seglist = runSegments( edgelist, minfrac, nonEdgeThresh, overwrite )
%
%  Run watershed segmentation on set of edge images
%  edgelist:  names of edges
%  nonEdgeThresh:  non-edge threshold, default 0.95
%                  will have no edges for density > nonEdgeThresh
%  minFrac:        can connect segments to larger segments if their shared
%                  boundary is at least this fraction of the smaller
%                  segment circumference, default 0.1
%  seglist:   list of segmentation images
%
%  Daniel Morris, Oct 2017
%
%  See also: scoreEdgeNetOnData, runEdgeSegNet, runNetEdges
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function seglist = runSegments( edgelist, minfrac, nonEdgeThresh, overwrite )

if nargin < 2, minfrac = 0.1; end
if nargin < 3, nonEdgeThresh = 0.98; end
if nargin < 4, overwrite = true; end

nim = numel(edgelist);
seglist = cell(1,nim);
for its=1:nim
    seglist{its} = strrep( edgelist{its}, '-edge','-segments');
end

fprintf('\n==> Running segmentation on %d edge images\n',  nim );

for its=1:nim
    %skip files that have already been calculated
    if ~overwrite && exist(seglist{its},'file')
        fprintf('Exists so skipping: %s\n', seglist{its} );
        continue;
    end
    
    fprintf('Reading %d of %d edge images, ', its, nim);
    edge = im2single( imread(edgelist{its}) );
    
    %Do watershed segment estimation using density image (from edge image):
    segLabels = edges2Segs( 1 - edge, nonEdgeThresh, minfrac );

    imwrite( segLabels2Colors( segLabels, true ), seglist{its} );
    fprintf('Saved segments\n');
end

end

