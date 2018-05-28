%imCmap
%
%  Class to easily annotate a grayscale image with colored pixels.  It
%  stores the image and so provides an easy way to use colormapBand.m
%  class.
%
%  Constructor:
%  imc = imCmap( img, nBands, colorFraction, colorbands)
%          img: image.  If color converts to grayscale.
%          nBands: number of color bands. 1 means only grayscale
%          colorFraction: Closer to 1 makes colors stronger,
%                         default: 0.25
%          colorbands:    [nBands-1 x 3] colors of bands
%                         default: jet(nBands-1)
%  imc.reset( img ) resets with new image (keeps same colormapping)
%  imc.reset( img, nBands, colorFraction ) resets with new image etc.
%  imc.setPixBand( pixinds, bandind)
%          adjusts the band for a subset of the pixels
%          pixinds: sets pixels img(pixinds)
%          bandind: sets them to this band (between 1 and nBands)
%  imc.show: displays with colorbands
%
%  Example usage:
%
%  img = imread('autumn.tif');
%  [u,v] = meshgrid((1:size(img,2))-size(img,2)/2, (1:size(img,1))-size(img,1)/2 );
%  pixSelect1 = u > v;
%  pixSelect2 = u+v > 0;
%  imc = imCmap( img, 4 );
%  imc.setPixBand( pixSelect1 & ~pixSelect2, 2 );  %band 2
%  imc.setPixBand( ~pixSelect1 & pixSelect2, 3 );  %band 3
%  imc.setPixBand( pixSelect1 & pixSelect2, 4 );  %band 4
%  imc.show
%
%  Daniel Morris, Dec 2015
%
%  See also: colormapBands
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
classdef imCmap < handle

    properties
        cimg      %image
        cmapBands %colormapBands class
    end
    
    methods
        
        function imc = imCmap( varargin )
            imc.cmapBands = colormapBands;
            if nargin>0
                imc.reset( varargin{:} );
            end
        end
        
        %varargin: colorFraction, colorbands
        function reset( imc, img, nBands, varargin  )
            imc.cmapBands.init( nBands, varargin{:} );
            imc.cimg = imc.cmapBands.im2Band( img );
        end
        
        function setPixBand( imc, pixinds, bandind )
            imc.cimg = imc.cmapBands.setPixBand( imc.cimg, pixinds, bandind );
        end
        
        function show( imc, varargin )
            imc.cmapBands.show( imc.cimg, varargin{:} );
        end
        
    end
    
end
