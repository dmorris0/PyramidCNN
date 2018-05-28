%segPixFeats
%
%  spf = segPixFeats( seglabels, pixopts, edgeDist, imEdgeSeg )
%
%  Calculate segment pixel features
%  seglabels: [M x N] image portion with pixel segment labels in range 1:P,
%                     with 0 being background
%  pixopts: 1xQ char array with the following chars:
%           'l': input edge label (i.e. seglabels)
%           'e': edge image (1 edges, 0 elsewhere)
%           'i': is-seg and unknown: 1 in segs, 0 elsewhere
%           'j': is-seg image: 1 in segs, 2 elsewhere
%           'n': is-seg and not edge (i.e. all segs > 0 but excluding edge
%                pixels
%           'b': 'e' pixels label 1, 'n' pixels label 2 (default)
%           'c': same as 'b', but background are label 3 (i.e clutter
%                rather than ignore)
%           'd': internal Euclidean distance from edge for all segments
%           '': then just returns spf = seglabels
%  edgeDist: max distance between pixels across boundary for them to be
%            labeled as edge pixels.  Default: 1, implies adjacent pixels
%            are edges but not diagonals.  If 1.5, then diagonals are also
%            edge pixels
%  imEdgeSeg: false (default) image boundaries are not edges
%             true, then segments that end at edge of the image have a
%             boundary along the image edge (but zero-seg does not)
%  spf: [M x N x Q] with channels being defined by pixopts
%
%  Daniel Morris, Oct 2017
%
%  See also: syntheticShapes, getImPatch, segs2edgeInside,
%            calcConPixPatchFilt, segEdgePR
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function spf = segPixFeats( seglabels, pixopts, edgeDist, imEdgeSeg )

if nargin < 4, imEdgeSeg = false; end
if nargin < 3, edgeDist = 1; end

pvals = unique( seglabels );

[M,N,P] = size(seglabels);
if P>1, error('seglabels must be M x N x 1'); end

if imEdgeSeg
    %then expand image with a zeros-boundary    
    nseg = zeros(M+2,N+2,'single');
    nseg(2:end-1,2:end-1) = seglabels;
    indist = zeros( M+2, N+2, 'single');
    outdist = zeros( M+2, N+2, 'single');
else
    nseg = seglabels;
    indist = zeros( M, N, 'single');
    outdist = zeros( M, N, 'single');
end

for its=1:numel(pvals)
    inpoly = (nseg == pvals(its));
    if pvals(its)==0
        outdist = bwdist(~inpoly);
    else
        indist = indist + bwdist(~inpoly);
    end
end

edge = indist + outdist <= edgeDist;
if imEdgeSeg
    indist = indist(2:end-1,2:end-1);
    edge = edge(2:end-1,2:end-1);
end

spf = zeros(M, N, numel(pixopts),'single');
for its=1:numel(pixopts)
    switch pixopts(its)
        case 'l'
            spf(:,:,its) = seglabels;
        case 'e'
            spf(:,:,its) = edge;
        case 'i'
            spf(:,:,its) = (seglabels > 0);
        case 'j'
            spf(:,:,its) = (seglabels > 0) + 2*(seglabels == 0);
        case 'n'
            spf(:,:,its) = (seglabels > 0) & ~edge;
        case 'b'
            spf(:,:,its) = edge + 2*((seglabels > 0) & ~edge);
        case 'c'
            layer = edge + 2*((seglabels > 0) & ~edge);
            layer(layer==0) = 3;
            spf(:,:,its) = layer;
        case 'd'
            spf(:,:,its) = indist;
        otherwise
            error(['Invalid pixopt: ' pixopts(its)]);
    end
end

end

