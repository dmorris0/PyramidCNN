%pixEdges
%
%  Class for defining and plotting edges between pixels.  An image has MxN
%  pixels and (M+1) x (N+1) edges.  These are indexed first by (M+1)*N
%  top-bottom edges, and then by M*(N+1) left-right edges.
%
%  Constructor:
%  pe = pixEdges( M, N )
%
%  Methods:
%  pe.init( M, N )
%  pe.plotSegBoundaries( seglabel, lineOrDot, varargin)
%                        segLabel: binary image
%                        lineOrDot: 'line' (default) 'dot' alternate plot
%                        It will plot the pixel boundaries.
%                        varargin: this is the line properties in a cell array ex: {'r-','linewidth',0.5};
%
%  einds = pe.getEdgeIndices( pixinds ): get edge indices given pixel
%               indices
%  pe.plotEdges( einds, sopt): plots edges
%  
%  Find boundary recall:
%   [BRecall, goodpnts, pMidGT, pMidTest] = pe.calcSegBoundaryRecall(... 
%                   segLabelGT, segLabelTest, edist, removeImageEdges)
%
%  Daniel Morris, Feb 2016
%
%  See also: plotSegBoundaries
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
classdef pixEdges < handle
   
    properties
        M  %number of rows of pixels in image
        N  %number of columns of pixels in image
    end
    
    methods
        %Constructor
        function pe = pixEdges(M, N)
            if nargin>0
                pe.init(M,N);
            end
        end
        
        function init(pe, M,N)
            pe.M = M;
            pe.N = N;
        end
        
        function pe = copy(ope)
            pe.M = ope.M;
            pe.N = ope.N;
        end
        
        %get edge indices for pixel(s) given its index
        %pixind: [m x n] pixel index
        %side:   'all' (default) returns all edges in [(m*n)x4] format [top bottom left right] 
        %        'left', 'right', 'top','bottom': returns just those edges
        %        in [m x n] matrix
        %einds: edges of pixels
        function einds = getEdgeIndices( pe, pixind, side )
            if nargin<3, side = 'all'; end
            if ~ischar(side), error('side must be string option'); end
            [ind,jnd] = ind2sub([pe.M pe.N], pixind );
            switch side
                case 'all'
                    top = sub2ind([pe.M+1, pe.N], ind, jnd);
                    bottom = sub2ind([pe.M+1, pe.N], ind+1, jnd);
                    left = sub2ind([pe.M, pe.N+1], ind, jnd) + (pe.M+1)*pe.N;
                    right = sub2ind([pe.M, pe.N+1], ind, jnd+1) + (pe.M+1)*pe.N;
                    einds = [top(:) bottom(:) left(:) right(:)];
                case 'left'
                    einds = sub2ind([pe.M, pe.N+1], ind, jnd) + (pe.M+1)*pe.N;
                case 'right'
                    einds = sub2ind([pe.M, pe.N+1], ind, jnd+1) + (pe.M+1)*pe.N;
                case 'top'
                    einds = sub2ind([pe.M+1, pe.N], ind, jnd);
                case 'bottom'
                    einds = sub2ind([pe.M+1, pe.N], ind+1, jnd);
                otherwise
                    error(['Invalid side: ' side]);
            end
        end
        
        %get edge endpoints coordinates
        %i.e. where an edge is the side of a square centered at a pixel
        %einds: edge indices, see pe.getEdgeIndices() above
        %p1: [2xNe] the start coordinates of each straight-line edge 
        %p2: [2xNe] the start coordinates of each straight-line edge
        %    p1 and p2 are separated by 1 pixel length
        function [p1, p2] = edgeEndCoords( pe, einds )
            tb = einds <= (pe.M+1)*pe.N; %top-bottom edges
            rl = ~tb; %right-left edges
            %Now get coordinates:
            [itb,jtb] = ind2sub([pe.M+1,pe.N], einds(tb));
            [irl,jrl] = ind2sub([pe.M,pe.N+1], einds(rl) - (pe.M+1)*pe.N);
            x1tb = jtb(:)'-0.5;
            x2tb = jtb(:)'+0.5;
            y1tb = itb(:)'-0.5;
            y2tb = itb(:)'-0.5;
            x1rl = jrl(:)'-0.5;
            x2rl = jrl(:)'-0.5;
            y1rl = irl(:)'-0.5;
            y2rl = irl(:)'+0.5;
            p1 = [x1tb x1rl; y1tb y1rl];
            p2 = [x2tb x2rl; y2tb y2rl];
        end

        %p: [2xNe] the middle coordinates of each straight-line edge 
        function p = edgeMiddleCoords( pe, einds )
            [p1,p2] = pe.edgeEndCoords( einds );
            p = (p1 + p2)/2;
        end
        
        %plot edges by index, either as lines or dots
        %einds: edge indices
        %lopt: 'line' (default) plots edges as lines
        %      'dot'  plots edges as dots in the middle of the edge
        %varargin: can optionally pass sopt cell into plotEdges function
        function plotEdgesByIndex(pe, einds, lopt, varargin)
            if nargin<3, lopt='line'; end
            %get edge end coordinates:
            switch lopt
                case 'line'
                    %add nans so that we can plot all edges together:
                    [p1,p2] = pe.edgeEndCoords( einds );
                    pe.plotEdges( p1, p2, varargin{:});
                case 'dot'
                    p = pe.edgeMiddleCoords( einds );
                    pe.plotEdges( p, [], varargin{:});
                otherwise
                    error('Invalid lopt');
            end
        end
        
        %p1: [2xNe] the start coordinates of each straight-line edge 
        %p2: [2xNe] the start coordinates of each straight-line edge
        %OR:
        %p1: [2xNe] the mid coordinates of each straight-line edge 
        %p2: []
        %sopt: cell array of line options
        function h = plotEdges(~, p1, p2, sopt)
            if nargin<3, p2 = []; end
            if nargin<4 || isempty(sopt)
                if ~isempty(p2)
                    sopt={'r-','linewidth',0.5};
                else
                    sopt={'r.','markersize',2};
                end
            end
            if ischar(sopt), sopt = {sopt}; end
            if ~isempty(p1) && size(p1,1)~=2, error('points should be 2xN');end
            if ~isempty(p2) && size(p2,1)~=2, error('points should be 2xN');end
            if isempty(p2)
                xn = p1(1,:);
                yn = p1(2,:);
            else
                xn = [p1(1,:);p2(1,:);nan(1,size(p1,2))];
                yn = [p1(2,:);p2(2,:);nan(1,size(p1,2))];
            end
            ah = ishold;
            if ~ah, hold on; end
            h = plot(xn(:),yn(:),sopt{:});
            if ~ah, hold off; end
            if nargout==0, clear h; end
        end
        
        
        %finds indices to edges (used by pe.plotEdges)
        %seglabel: [MxN] image with label for each pixel giving its segment
        function einds = boundaryEdgeIndices(pe, seglabel, includeImageBoundary)
            if isempty(pe.M), error('pixEdges not initialized'); end
            if size(seglabel,1)~= pe.M || size(seglabel,2)~=pe.N
                error('seglabel should be size M x N');
            end
            if nargin < 3, includeImageBoundary = true; end
            tb = false(pe.M+1,pe.N);
            tb(2:end-1,:) = seglabel(2:end,:) ~= seglabel(1:end-1,:) & ~isnan( seglabel(2:end,:) + seglabel(1:end-1,:) );
            if includeImageBoundary
                tb([1 end],:) = true;
            end
            lr = false(pe.M,pe.N+1);
            lr(:,2:end-1) = seglabel(:,2:end) ~= seglabel(:,1:end-1) & ~isnan( seglabel(:,2:end) + seglabel(:,1:end-1) );
            if includeImageBoundary
                lr(:,[1 end]) = true;
            end
            einds = [find(tb(:)') find(lr(:)')+(pe.M+1)*pe.N];
        end
        
        %in order to compare my method with other methods, I need to be
        %able to plot their results in the same way.  This function takes
        %in a segmentation label image (assumed to be size MxN) with a
        %label for each pixel.  It will plot the pixel boundaries.
        %varargin: this is the line properties in a cell array ex: {'r-','linewidth',0.5};
        function plotSegBoundaries(pe, seglabel, lineOrDot, sopt, includeImageBoundary )
            if nargin < 4, sopt = {}; end
            if nargin < 5, includeImageBoundary = true; end
            if isempty(pe.M) || pe.M ~= size(seglabel,1) || pe.N ~= size(seglabel,2)
                pe.init( size(seglabel,1), size(seglabel,2) );
            end
            if nargin<3, lineOrDot = 'line'; end
            pe.plotEdgesByIndex( pe.boundaryEdgeIndices(seglabel, includeImageBoundary), lineOrDot, sopt );
        end
        
        %Finds boundary recall, namely fraction of edges in ground truth
        %using centers of each edge point
        %that are within edist of an edge in the test segments
        %segLabelGT:   [MxN] labels for segments
        %segLabelTest: [MxN] labels for testing segments
        %edist:        scalar: euclidean distance in pixels
        %BRecall:       Boundary Points Recall  fraction of ground truth points that are withing
        %                  edist of any test point
        %goodpnts: [Ng x 1] logical array identifying which ground truth
        %                   points are good
        %pMidGT:   [2 x Ng] edge midpoints of ground truth
        %pMidTest: [2 x Nt] edge midpoints of test segments
        %removeImageEdges: true (default) removes boundary pixels on the
        %                  edge of the image as these should not be counted
        %                  in the boundary recall
        function [BRecall, goodpnts, pMidGT, pMidTest] = calcSegBoundaryRecall(pe, segLabelGT, segLabelTest, edist, removeImageEdges)
            if nargin < 5, removeImageEdges = true; end
            pMidGT = pe.edgeMiddleCoords( pe.boundaryEdgeIndices( segLabelGT ) ); 
            pMidTest = pe.edgeMiddleCoords( pe.boundaryEdgeIndices( segLabelTest ) ); 
            if removeImageEdges
                %now remove image boundary pixels:
                pMidGT = pe.removeBoundaryEdgePix( pMidGT );
                pMidTest = pe.removeBoundaryEdgePix( pMidTest );
            end
            %Finally calculate the Boundary Recall:
            [BRecall,goodpnts] = pe.calcBRecall(  pMidGT, pMidTest, edist);
        end
        
        %find Precision-Recall curve for an edge image
        %Note: this is not too useful since it calls calcBRecall too may
        %times.  Instead use measure segEdgePR
        function [rlist, plist, ap] = calcEdgeImPRcurve(pe, segLabelGT, edgeIm, edist, threshlist )
            if nargin < 5
                threshlist = linspace(min(edgeIm(:)), max(edgeIm(:)), 100 );
            end
            %get target edge points which we want to estimate:
            pMidGT = pe.edgeMiddleCoords( pe.boundaryEdgeIndices( segLabelGT ) ); 
            pMidGT = pe.removeBoundaryEdgePix( pMidGT );
            [u,v] = meshgrid(1:pe.N, 1:pe.M);
            %which portion of image is within edist of edge
            validRegions = segPixFeats( segLabelGT+1, 'd' ) < (edist + 1);  %add 1 to edist since this is actually distance from pixels
            [rlist,plist] = arrayfun( @(x) calcPR( x ), threshlist );
            
            %calculate average precision by summing area under curve:
            barheight = 0.5*(plist(1:end-1)+plist(2:end));
            barwidth = rlist(2:end) - rlist(1:end-1);
            ap = barheight * barwidth';
            
            function [r,p] = calcPR( t )
                egood = edgeIm(:)' > t & validRegions(:)';
                ecoords = [u( egood ); v( egood )];
                %find recall of boundary points:
                [r, gb, ge] = pe.calcBRecall( pMidGT, ecoords, edist );
                %find bad edge points that are in non-background segs:
                ebad = edgeIm > t & segLabelGT > 0;
                etaken = find(egood);
                etaken = etaken(ge);
                ebad(etaken) = false;  %exclude points that are close enough to edges
                fp = sum( ebad(:) );
                %Precision = TP / (TP + FP)
                tp = sum(gb);
                p = tp / (tp + fp);
                %This measure is a little strange as fp refers to edge
                %pixels and tp refers to mid-points between pixels on
                %boundary.  Also all pixels within edist of a tp are
                %removed from detections.  I think that is correct to do. 
            end
        end
            
            
        
        %see calcSegBoundaryRecall for args
        %mask: save size as segLabelGT and segLabelList with values > 0 t
        %col = [0 1 0;    %seg GT good
        %       1 0 0;    %seg GT bad
        %       1 1 0;    %test good
        %       1 0 1];   %test bad
        %       nan: then no plotting
        function [BRecall, BPrecision, nGT, nTest] = plotBoundaryWithGT( pe, segLabelGT, segLabelTest, edist, mask, col )
            if nargin < 5, mask = []; end
            if nargin < 6 || isempty(col)
                col = [0 1 0; %seg GT good
                    1 0 0;    %seg GT bad
                    1 1 0;    %test good
                    1 0 1;   %test bad inside
                    1 1 1];   %test bad outside
            elseif isscalar(col) && isnan(col)
                col = nan(4,3);
            else
                if size(col,1) < 4 || size(col,2) < 3
                    error('col must be 4 x 3');
                end
            end
            if ~isempty(mask)
                segLabelGT(mask==0) = nan; 
                segLabelTest(mask==0) = nan;
            end
            pMidGT = pe.edgeMiddleCoords( pe.boundaryEdgeIndices( segLabelGT ) ); 
            pMidTest = pe.edgeMiddleCoords( pe.boundaryEdgeIndices( segLabelTest ) ); 
            pMidGT = pe.removeBoundaryEdgePix( pMidGT );
            pMidTest = pe.removeBoundaryEdgePix( pMidTest );
            
            %Finally calculate the Boundary Recall:
            [BRecall,goodpnts,goodtestpnts, BPrecision] = pe.calcBRecall(  pMidGT, pMidTest, edist);
            nGT = numel(goodpnts);
            nTest = numel(goodtestpnts);
            badTest = find(~goodtestpnts);
            badTestInds = sub2ind([pe.M pe.N], floor(pMidTest(2,~goodtestpnts)), floor(pMidTest(1,~goodtestpnts)) );
            badInside = segLabelGT(badTestInds) > 0;
            ah = ishold;
            if ~ah, hold on; end
            if ~isnan(col(1,1)), plot(pMidGT(1,goodpnts), pMidGT(2,goodpnts), '.', 'color', col(1,:) ); end
            if ~isnan(col(2,1)), plot(pMidGT(1,~goodpnts), pMidGT(2,~goodpnts), '.', 'color', col(2,:) ); end
            if ~isnan(col(3,1)), plot(pMidTest(1,goodtestpnts), pMidTest(2,goodtestpnts), '.', 'color', col(3,:) ); end
            %if ~isnan(col(4,1)), plot(pMidTest(1,~goodtestpnts), pMidTest(2,~goodtestpnts), '.', 'color', col(4,:) ); end
            if ~isnan(col(4,1)), plot(pMidTest(1,badTest(badInside)), pMidTest(2,badTest(badInside)), '.', 'color', col(4,:) ); end            
            if ~isnan(col(5,1)), plot(pMidTest(1,badTest(~badInside)), pMidTest(2,badTest(~badInside)), '.', 'color', col(5,:) ); end
            if ~ah, hold off; end
        end
        
        %output for each image pixel the distance to the closest boundary
        %segLabels: M x N segment labels with boundaries defined between
        %them
        function distim = calcDistToClosestBoundary(pe, segLabels, removeImageEdges)
            if nargin < 3, removeImageEdges = true; end
            pMidGT = pe.edgeMiddleCoords( pe.boundaryEdgeIndices( segLabels ) ); 
            if removeImageEdges
                %now remove image boundary pixels:
                pMidGT = pe.removeBoundaryEdgePix( pMidGT );
            end
            [u,v] = meshgrid(1:size(segLabels,2),1:size(segLabels,1));
            distim = reshape( min( pdist2( [u(:) v(:)], pMidGT' ), [], 2), size(segLabels,1), size(segLabels,2));
        end
        
        %remove pixels that are on the boundary of the image as these
        %should not be included in the statistics
        function pts = removeBoundaryEdgePix( pe, pts )
            bpix = pts(1,:) == pe.N+0.5 | pts(1,:) == 0.5 | pts(2,:) == 0.5 | pts(2,:) == pe.M+0.5;
            pts(:,bpix) = [];
        end
        
        %Calculates boundary recall
        %gtpnts:   [2 x Ng] ground truth points, each row is [x y]
        %testpnts: [2 x Nt] test points with x,y rows
        %BRecall:  scalar, fraction of ground truth points that are withing
        %                  edist of any test point
        %goodpnts: [Ng x 1] logical array identifying which ground truth
        %                   points are good
        function [BRecall, goodpnts, goodtestpnts, BPrecision] = calcBRecall(~, gtpnts, testpnts, edist )
            if size(gtpnts,1) > 2 || size(testpnts,1) > 2
                error('gtpnts and testpnts must be 2 x Ng, Nt');
            end
            pd = pdist2( gtpnts', testpnts' );
            goodpnts = min( pd, [], 2) <= edist + 1e-5; %add a delta to avoid roundoff error when passing in integer distances
            BRecall = sum(goodpnts) / numel(goodpnts);
            if nargout > 2
                goodtestpnts = min( pd, [], 1)' <= edist + 1e-5;
                BPrecision = sum(goodtestpnts) / numel(goodtestpnts);
            end
        end
        
        %plots edges as dots for ground truth segments in:
        %red if within edist of a superpixel
        %green if not
        %plots segLabelTest edges in magenta 
        function BRecall = plotSegBoundaryRecall(pe, segLabelGT, segLabelTest, edist )
            ah = ishold;
            if ~ah, hold on; end
            [BRecall, goodpnts, pMidGT] = calcSegBoundaryRecall(pe, segLabelGT, segLabelTest, edist);
            %Indices of gt edges:
            eIndsGT = pe.boundaryEdgeIndices( segLabelGT );
            [pEnd1,pEnd2] = pe.edgeEndCoords( eIndsGT );
            %now plot good and bad in different colors.  Plot both dots in
            %the middle of edges and line for full edges:
            pe.plotEdges( pMidGT(:,goodpnts), [], {'r.'} );
            pe.plotEdges( pEnd1(:,goodpnts), pEnd2(:,goodpnts), {'r:'} );
            pe.plotEdges( pMidGT(:,~goodpnts), [], {'g.'} );
            pe.plotEdges( pEnd1(:,~goodpnts), pEnd2(:,~goodpnts), {'g:'} );
            %Finally plot segLabelTest
            pe.plotSegBoundaries( segLabelTest, 'line', {'b-','color',[.5 0 1]} );
            if ~ah, hold off; end;
            title(['Boundary Recall: ' num2str(BRecall,'%.3f')]);
        end
        
    end
        
end
