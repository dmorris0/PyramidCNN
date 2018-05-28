%plotSegBoundaries
%
%  plotSegBoundaries( seglabels )
%  plotSegBoundaries( seglabels, lopt )
%  plotSegBoundaries( seglabels, lopt, includeImageBoundary)
%
%  Plots segment boundaries on current image
%  segLabels: labeled image with one label per pixel
%  lopt: line properties in a cell array ex: {'k-','linewidth',2}
%
%  Daniel Morris, Oct 2016
%
%  See also: pixEdges, plotSegs, plotColSegBoundaries
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function plotSegBoundaries( seglabels, varargin )

pe = pixEdges;
pe.plotSegBoundaries( seglabels, 'line', varargin{:} );

end

