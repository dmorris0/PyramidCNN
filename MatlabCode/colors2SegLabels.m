%colors2SegLabels
%
%  segLabels = colors2SegLabels( colorSegs )
%
%  Convert colored segment image into numbered label image with a different
%  number for each segment in range 1 to N-segs.
%  colorSegs: [M x N x P] must be uint8
%  segLabels: [M x N] numbered segments
%
%  Daniel Morris, Oct 2016
%
%  See also: segLabels2Colors
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function segLabels = colors2SegLabels( colorSegs, doSort )

if nargin <2, doSort = true; end

if ~isa(colorSegs,'uint8'), error('Only operates on uint8 images'); end

[M,N,P] = size(colorSegs);

%first get single number for each color:
labels = zeros(M,N);
for its = 1:P
    labels = labels*256 + double(colorSegs(:,:,its));
end

if doSort
    L = unique(labels)';
    L(L==0) = [];  %zero colors remain zero
    inds = 1:numel(L);
    
    segLabels = zeros(M,N);
    arrayfun(@(x,y)setLabelIndex(x,y), L, inds );
else
    segLabels = labels;
end

    function setLabelIndex( lval, ind )
        segLabels(labels==lval) = ind;
    end


end
