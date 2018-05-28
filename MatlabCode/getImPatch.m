%getImPatch
%
%  imp = getImPatch( img, bb )
%  imp = getImPatch( img, xymid, xywid )
%      
%  Returns image patch.  Fills any region outside the image boundary with
%  zeros
%  img: [M x N x P] image
%  bb:  [xmin xmax ymin ymax]
%  xymid: [xm ym] center of rectangle.  If widths are even, then this
%          should be on the half-pixel to be properly centered
%  xywid: [xw yw] widths of rectangle (and output imp)
%  Note: all units rounded to closest integer pixel coords
%
%  Daniel Morris, Oct 2017
%
%  See also: im2grid, fixImDim
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function imp = getImPatch( img, xymid, xywid )

if nargin < 3
    bb = round( xymid );
    xywid = [bb(2)-bb(1)+1 bb(4)-bb(3)+1];
else
    xymin = round( xymid - (xywid-1)/2 );
    xywid = round( xywid );
    bb = round( [xymin(1)+[0 xywid(1)-1] xymin(2)+[0 xywid(2)-1]] ); 
end

imp = zeros(xywid(2),xywid(1),size(img,3),size(img,4),'like',img);

[M,N,~] = size(img);
xi = bb(1):bb(2);
xin = xi > 0 & xi <= N;
yi = bb(3):bb(4);
yin = yi > 0 & yi <= M;
imp(yin,xin,:) = img(yi(yin),xi(xin),:);

end
