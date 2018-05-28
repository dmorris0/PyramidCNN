%fixImDim
%
%  fimg = fixImDim( img, factor, doExpand )
%      
%  Fix image dimensions -- ensure that they are a multiple of factor.
%  img:    [M x N x P x Q]
%  factor: new dimensions, nM and nN will be divisible by factor
%  fimg:   [nM x nN x P x Q]
%  doExpand: true (default): nM and nN >= M and N respectively
%                 zero-padding is used
%            false: nM and nN <= M, N respectively
%
%  Daniel Morris, Nov 2017
%
%  See also: im2grid, getImPatch
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function img = fixImDim( img, factor, doExpand )

if nargin < 3, doExpand = true; end
if factor < 1, error(['factor must be >=1, input: ' num2str(factor)]); end
csize = [size(img,2) size(img,1)];
if doExpand
    nsize = ceil( csize/factor ) * factor;
else
    nsize = floor( csize/factor ) * factor;
end    
if any(csize~=nsize)
    %then got to resize image
    img = getImPatch(img,[1 nsize(1) 1 nsize(2)]);
end

end
