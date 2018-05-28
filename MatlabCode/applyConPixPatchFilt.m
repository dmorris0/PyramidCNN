%applyConPixPatchFilt
%
%  offsegchan = applyConPixPatchFilt( segim, cppFilt )
%
%  Apply Connected Pixel Patch filter to an image.  See notes: Learning / 17-10-25
%  segim:    [M x N x 1 x P] class single image
%  cppFilt:  [n x n x nchan] An offset filter that for each channel will
%            shift the image a certain number of pixels within the radius
%            See: calcConPixPatchFilt.m for calculation of it
%  offsegchan: [M x N x nchan x P] 
%
%  Daniel Morris, Oct 2017
%
%  See also: calcConPixPatchFilt, getImPatch, segs2edgeInside,
%            segPixFeats, tryGPU
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function offsegchan = applyConPixPatchFilt( segim, cppFilt  )

if ~isa(segim,'single')
    error('Segim must be singles');
end
assert( size(segim,3) == 1 );

[sy,sx,sz] = size(cppFilt);

padval = floor(sx/2);

%shift dimensions of cppFilt since input is 1 channel (in dimension 3) and output is sz
%channels in dimension 3
if exist('vl_nnconv','file') ~= 3
    error('It is necessary to install MatConvNet (http://www.vlfeat.org/matconvnet/), build it (CPU version is fine) and add the vl_nnconv mex folder to the path');
end
offsegchan = vl_nnconv( segim, reshape( cppFilt, [sy, sx, 1, sz] ), [], 'pad', padval );

end

