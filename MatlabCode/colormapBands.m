%colormapBands
%
%  Class to break indexed image into color bands
%  Typically applies to grayscale image to enable coloring of certain
%  pixels while maintaining their rough intensity.
%  See: imCmap for a simple way to use this.
%
%  Constructor:
%  cb = colormapBands( nBands, colorFraction )
% 
%  Methods:
%  cb.init( nBands, colorFraction ):  initialize bands with default colors.  Can
%                      optionally vary these colors
%                      colorFraction: Closer to 1 makes colors stronger,
%                                     default: 0.25
%                      colorbands:    [nBands-1 x 3] colors of bands
%                                     default: jet(nBands-1)
%  cb.setColoredGrayscaleMap( colorFraction, colorbands ): can set
%                      non-default color bands
%  cimg = cb.im2Band( img, imgrange ): transforms image into a
%                      single-band range
%  cimg = cb.setPixBand( cimg, pixinds, bandind): adjusts the band for a
%                      subset of the pixels
%                      pixinds: sets pixels img(pixinds)
%                      bandind: sets them to this band (between 1 and
%                      nBands)
%  cb.show( cimg ):    displays image with colorbands
%                      cimg must be indexed into the colormap 
%  cb.demo( img ):  demo for using this class
%
%  Example usage:
%
%  cb = colormapBands( 4 );                          %want 4 color bands
%  img = imread('autumn.tif');
%  cimg = cb.im2Band( img );                         %All pixels initially in band 1
%  [u,v] = meshgrid((1:size(img,2))-size(img,2)/2, (1:size(img,1))-size(img,1)/2 );
%  pixSelect1 = u > v;
%  pixSelect2 = u+v > 0;
%  cimg = cb.setPixBand( cimg, pixSelect1 & ~pixSelect2, 2 );  %band 2
%  cimg = cb.setPixBand( cimg, ~pixSelect1 & pixSelect2, 3 );  %band 3
%  cimg = cb.setPixBand( cimg, pixSelect1 & pixSelect2, 4 );   %band 4
%  cb.show( cimg );
%
%  Daniel Morris, Sep 2015
%
%  See also: imCmap
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
classdef colormapBands < handle

    properties
        nBands      %number of bands
        nPixPerBand %number of levels per band
        nCol        %total number of levels (nBands*nPixPerBand)
        cmap        %[nCol x 3] colormap
        cFraction
        colorbands
    end
    
    methods
        
        function cb = colormapBands( varargin )
            cb.cFraction = 0.25;
            cb.colorbands = [];
            if nargin<1
                cb.init(1);
            else
                cb.init( varargin{:}); 
            end
        end
        
        function init(cb, nBands, colorFraction, colorbands)
            if nargin<3, colorFraction = cb.cFraction; end
            if nargin<4, colorbands = cb.colorbands; end
            cb.nBands = nBands;
            cb.nPixPerBand = floor(256/nBands);
            cb.nCol = cb.nBands * cb.nPixPerBand;
            cb.setColoredGrayscaleMap( colorFraction, colorbands );
        end
      
        %create a new band image, and optionally display it
        %cimg = cb.create( img, pmask1, pmask2, ... )  
        %will re-initializes cb with 1 band for image 1 band per mask
        %Note; later masks override earlier ones
        function cimg = create(cb, img, varargin )
            cb.init( nargin-1 );
            if islogical(img)
                cimg = cb.im2Band( 255*double(img) );
            else
                cimg = cb.im2Band( mean( img, 3) );
            end
            for its=1:nargin-2
                if ~islogical(varargin{its}), error('pixel masks must be logicals'); end
                cimg = cb.setPixBand(cimg, varargin{its}, its+1);
            end
            if nargout==0
                cb.show( cimg );
                clear cimg
            end
        end

        %set grayscale map with colored bands
        %colorFraction: amount of color in bands, default 0.25 (between 0
        %               and 1)
        %               Closer to 1 makes colors stronger
        %colorbands:    [nBands-1 x 3] colors of bands
        %               default: jet(nBands-1)
        function setColoredGrayscaleMap( cb, colorFraction, colorbands )
            if nargin<2, colorFraction = cb.cFraction; end
            if nargin<3 || isempty(colorbands), colorbands = myDefaultColorMap(cb.nBands-1); end
            if size(colorbands,1) ~= cb.nBands-1
                error(['Should have ' num2str( cb.nBands-1 ) ' rows in colorbands']);
            end
            cb.cmap = zeros( cb.nBands * cb.nPixPerBand, 3);
            cb.cmap(1:cb.nPixPerBand,:) = gray(cb.nPixPerBand);
            if cb.nBands > 1
                bscale = 1-colorFraction;
                for its=2:cb.nBands
                    boffset = ones(cb.nPixPerBand,1) * colorFraction * colorbands(its-1,:);
                    cb.selectBandColors( its, gray(cb.nPixPerBand)*bscale + boffset );
                end
            end
            function dmap = myDefaultColorMap( n )
                %typically few number of bands so select furthest apart
                %ones first:
                if n<=4
                    dmap = [1 .4 0;.2 .3 1;.2 1 0;0 .5 .6];
                    dmap = dmap(1:n,:);
                else
                    dmap = jet(n);
                end
            end
        end
        
        function selectBandColors( cb, ind, bcmap )
            if size(bcmap,1) ~= cb.nPixPerBand
                error(['Band colormap must have ' num2str(cb.nPixPerBand) ' rows']);
            end
            if size(bcmap,2) ~=3, error('colormaps must be nx3'); end;
            if max(bcmap(:))>1 || min(bcmap(:)) < 0
                error('Colormaps must be in range 0 to 1');
            end
            indlist = (ind-1)*cb.nPixPerBand+1:ind*cb.nPixPerBand;
            cb.cmap(indlist,:) = bcmap;
        end
    
        %transform image img so that its values are integers within the
        %range of 1 band
        %imgrange: [minval maxval], default: [min(img(:)) max(img(:))]
        function cimg = im2Band( cb, img, imgrange )
            if size(img,3)>1
                img = rgb2gray(img);  %make grayscale
            end
            %make sure img is a double:
            img = double(img);
            if nargin<3 || isempty(imgrange)
                imgrange = [min(img(~isnan(img) & ~isinf(img))) max(img(~isnan(img) & ~isinf(img)))];
            end
            if sum(isnan(imgrange))>0
                warning('imgrange has nans');
            end
            cimg = (img-imgrange(1)) / (imgrange(2)-imgrange(1));
            cimg(cimg<0) = 0;
            cimg(cimg>1) = 1;
            cimg = round(cimg*(cb.nPixPerBand-1)+1);
        end

        %offset pixels to diffent band
        function cimg = setPixBand( cb, cimg, pixinds, bandind)
            if bandind<1 || bandind>cb.nBands, error('bandind outside range'); end;
            cimg(pixinds) =  cb.nPixPerBand*(bandind-1) + mod(cimg(pixinds)-1,cb.nPixPerBand) + 1;
        end
        
        %displays image with colormap
        function show( cb, cimg, useImfun )
            if nargin<3, useImfun = true; end
            if useImfun
                im(cimg,[1 cb.nCol]) %necessary to add range otherwise im rescales the colormap
            else
                imagesc(cimg,[1 cb.nCol]) %necessary to add range otherwise imagesc rescales the colormap
            end
            colorbar('off');
            colormap(cb.cmap)
            title('');
        end            
        
        function demo( cb, img )
            cb.init( 4 ); %use 4 color bands
            if nargin<2, img = imread('autumn.tif'); end;
            cimg = cb.im2Band( img );                         %all pixels initially in band 1
            [u,v] = meshgrid((1:size(img,2))-size(img,2)/2, (1:size(img,1))-size(img,1)/2 );
            pixSelect1 = u > v;
            pixSelect2 = u+v > 0;
            cimg = cb.setPixBand( cimg, pixSelect1 & ~pixSelect2, 2 );  %band 2
            cimg = cb.setPixBand( cimg, ~pixSelect1 & pixSelect2, 3 );  %band 3
            cimg = cb.setPixBand( cimg, pixSelect1 & pixSelect2, 4 );   %band 4
            cb.show( cimg );
        end
    end
end

