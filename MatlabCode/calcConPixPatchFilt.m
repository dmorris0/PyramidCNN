%calcConPixPatchFilt
%
%  [cppFilt, xyoff] = calcConPixPatchFilt( radius, includeZeroOffset, onlyHalf )
%  [cppFilt, xyoff] = calcConPixPatchFilt( xyoff )
%
%  Calculate Connected Pixel Patch filter.  See notes: Learning / 17-10-25
%  radius:   radius of disk
%  includeZeroOffset: false (default) excludes zero-offset channel in cppFilt
%                     true: then first channel is zero-offset
%  onlyHalf: false (default): the full disk around each pixel
%            true: then only portion of disk where pixels have lower index
%  cppFilt:  [n x n x nchan] An offset filter that for each channel will
%            shift the image a certain number of pixels within the radius
%  xyoff:    [nchan x 2] the xy offset obtained on a target image when
%            applying each filter, see applyConPixPatchFilt
%  Note: the channels are sorted so that magnitude of offset increases
%  (except when xyoff is input)
%
%  Daniel Morris, Oct 2017
%
%  See also: applyConPixPatchFilt, getImPatch, segs2edgeInside,
%            segPixFeats, seg2ConPixPatch, 
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function [cppFilt, xyoff] = calcConPixPatchFilt( radius, includeZeroOffset, onlyHalf )

if nargin < 2, includeZeroOffset = false; end
if nargin < 3, onlyHalf = false; end

if isscalar(radius)
    if radius > 1
        hn = floor(radius-eps(radius));  %for larger patches (with integer radii) it does not make sense for the outside row/column to have only one 1 in it
    else
        hn = floor(radius);
    end
    n = 2*hn+1;
    [x,y] = meshgrid((1:n)-(hn+1), (1:n)-(hn+1) );
    mask = x.^2 + y.^2 <= radius^2;
    if ~includeZeroOffset
        %if we don't want zero offset, zero this:
        mask(hn+1,hn+1) = false;
    end
    if onlyHalf
        %mark portion with indices greater than center as false:
        ind = (x-1)*n+y;
        mask( ind > ind(hn+1,hn+1) ) = false;
    end
    nchan = sum(mask(:));
    xyoff = zeros(nchan,2);
    inc = 0;
    for its = 1:n
        for jts = 1:n
            if mask(jts,its)
                inc = inc+1;
                xyoff(inc,:) = (hn+1)-[its, jts];
            end
        end
    end
    
    %now let's sort by offset magnitude, with smallest first:
    r = sum( xyoff.^2, 2);
    [~,sind] = sort(r);
    xyoff = xyoff(sind,:);
else
    xyoff = radius;
    assert(size(xyoff,2)==2);
end

cppFilt = getFilt( xyoff );

end

function cppFilt = getFilt( xyoff )

hn = max(abs(xyoff(:)));
n = 2*hn+1;
nchan = size(xyoff,1);
cppFilt = zeros(n,n,nchan,'single');
for its=1:nchan
    cppFilt(hn+1-xyoff(its,2), hn+1-xyoff(its,1), its) = single(1);
end
end

%% Test
function test
%copy and paste the below to test calcConPixPatchFilt as well as
%halfNb2AllNb:
seg = ones(3,3,'single');

[FcppFilt, Fxyoff] = calcConPixPatchFilt( 1, false, false );
Fseg = applyConPixPatchFilt( seg, FcppFilt);
 
% Is this the same as the below:
[HcppFilt, Hxyoff] = calcConPixPatchFilt( 1, false, true );
Hseg = applyConPixPatchFilt( seg, HcppFilt);
[FHseg, FHxyoff] = halfNb2AllNb( Hseg, Hxyoff );

%To compare them one must sort them:
p1 = Fxyoff(:,1)*4 + Fxyoff(:,2);
p2 = FHxyoff(:,1)*4 + FHxyoff(:,2);

[~,s1] = sort(p1);
[~,s2] = sort(p2);

SFseg = Fseg(:,:,s1);
SFxyoff = Fxyoff(s1,:);

SFHseg = FHseg(:,:,s2);
SFHxyoff = FHxyoff(s2,:);

cdiff(SFseg,SFHseg)
cdiff(SFxyoff, SFHxyoff)


end
