%segLabels2Colors
%
%  cSegLabels = segLabels2Colors( segLabels )
%  cSegLabels = segLabels2Colors( segLabels, doRGB )
%  cSegLabels = segLabels2Colors( segLabels, img )
%
%  Takes a segment label image and outputs another with random labels in
%  range 1 to 255.  That way it can be displayed with a colormap.  This can have repeated colors. 
%  Can also convert to a rgb image with 3 layers. Note: this will not have
%  repeated colors unless there are more than 255^3 labels.
%  segLabels that are 0 remain zero, remainder are non-zero.
%  doRGB: false (default) converts to single layer output
%         true: converts to 3-layer color output
%  img:   [M x N x 3] in this case each label will be the average color of
%         the segment in img
%
%  Daniel Morris, Oct 2016
%
%  See also: colors2SegLabels
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function cSegLabels = segLabels2Colors( segLabels, doRGB )

if nargin<2, doRGB = false; end

if isscalar(doRGB)
    if doRGB
        %renumber in range 1 to 255^3
        segLabels = segLabels * floor( 255^3 / max(segLabels(:)) );
        cSegLabels = uint8( cat(3, segLabels2Colors( segLabels, false ), ...
            segLabels2Colors( floor(segLabels/255), false ), ...
            segLabels2Colors( floor(segLabels/255^2), false ) ) );
    else
        goodLabels = segLabels > 0;
        ng = sum(goodLabels(:));
        nSegLabels = mod(segLabels(goodLabels)-1,254)+1;  %put in range 1 to 255
        
        %create transformation:
        T = sparse(1:255,randperm(255),ones(1,255),255,255);
        
        M = sparse(nSegLabels(:)',1:ng,ones(1,ng), 255, ng );
        
        [tSegLabels,~] = find( T*M );
        
        cSegLabels = segLabels;
        cSegLabels(goodLabels) = tSegLabels;
    end
else
    img = doRGB;
    if size(img,1) ~= size(segLabels,1) || size(img,2) ~= size(segLabels,2)
        error('img must be same size as segLabels');
    end
    segs = unique(segLabels);
    outim = cell(1,size(img,3));
    for its=1:size(img,3)
        layer = img(:,:,its);
        outlayer = zeros(size(segLabels), 'like', img );
        arrayfun( @(x) setSeg(x), segs );        
        outim{its} = outlayer;
    end
    cSegLabels = cat(3, outim{:} );
end


    function setSeg( x )
    pix = segLabels == x;
    outlayer( pix ) = mean( layer( pix ) );
    end

end
    
