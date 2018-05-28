%pixNeighbor
%
%  Class to find indices of pixel neighbors.  Useful for building sparse
%  matrices of edges.  Assumes image pixels are converted into a vertical
%  vector [M*N x 1].  Finds list of indices to desired neighbor defined as
%  an offset: deltaRow, deltaCol.  Takes into account that pixels don't
%  have neighbors across boundaries of the image.
%
%  Constructor:
%   pn = pixNeighbor(M, N)
%
%   [plist, nlist] = pn.calcInds( deltaRows, deltaCols, wt):
%                               neighbor defined by deltaRow, deltaCol
%                               plist: indices of pixels with neighbors
%                               nlist: indices of corresponding neighbors
%   sa = pn.makeSparse( deltaRows, deltaCols, wt)
%                               deltaRows, deltaCols: [1xQ] list of
%                               neighbors for each pixel
%                               wt:  [1xQ] weight for each neighbor type, default: ones(1,Q)
%   rc = pn.makeEdges( deltaRows, deltaCols, wt)
%                               Similar to makeSparse, but instead of
%                               returing a sparse array, returns indices of
%                               connected edges in in row-column format
%                               rc: [rlist;clist] [2 x Nrc] 
%
%  Example usage (from pixMerger.m):
%            pn = pixNeighbor(pm.M, pm.N);
%            switch pm.params.connected,
%                case 8,
%                    pm.adjacent = pn.makeSparse([0 1 -1 1],[1 0 1 1],[1 1 sqrt(2)/2 sqrt(2)/2]); %8-connected with diagonals weighted sqrt(2)/2
%                case 4,
%                    pm.adjacent = pn.makeSparse([0 1],[1 0],[1 1]);  %4 connected
%                otherwise,
%                    error('Invalid connection, must be 4 or 8');
%            end
%
%
%  Daniel Morris, Oct 2015
%
%  See also: edgeNeighbor, pixMerger
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
classdef pixNeighbor

    properties
        M   %number of rows
        N   %number of columns
    end
    
    methods
        %constructor:
        function pn = pixNeighbor(M, N)
            pn.M = M;
            pn.N = N;
        end
        
        %Find indices of neighbors
        %dRow, dCol: row and column offset defining a neighbor
        %Note: can also work if dRow and dCol are vectors of offsets
        %dRows, dCols: [1xQ] row and column offset of Q neighbors
        %wt: [1xQ] weight for each neighbor type, default: ones(1,Q)
        %    OR {1 x Q} each element is [M x N], and wlist is found by
        %    using corresponding plist to look up values in this.  
        %plist: [1xP] indices of pixels having a neighbor
        %nlist: [1xP] indices of corresponding neighbors
        %wlist: [1xP] weights for each edge
        %pixMask: [] then all pixels in MxN
        %         [MxN] logical, then restricts pixels and neighbors to
        %         those both with pixMask(plist) == 1 and pixMask(nlist) == 1
        %Note: when using BK_SetNeighbors, dRow + M*dCol should be >=0 as
        %only the upper triangular  portion is considered
        function [plist, nlist, wlist] = calcInds( pn, dRow, dCol, wt, pixMask )
            if nargin<4, wt = ones(1,numel(dRow)); end;
            if nargin<5, pixMask = []; end;
            if numel(dRow)>1
                if numel(dRow) ~= numel(dCol) || numel(dRow) ~= numel(wt)
                    error('dRow, dCol and wt must be vectors of equal length');
                end
                [pc,nc,wc] = arrayfun( @(x,y,z) calcInds(pn, x, y, z, pixMask), dRow, dCol, wt, 'uniformoutput', false);
                plist = horzcat(pc{:});
                nlist = horzcat(nc{:});
                wlist = horzcat(wc{:});
                return;
            end
            %pix and nb will have the same size blocks of ones, just offset
            pix = ones(pn.M,pn.N);
            nb = ones(pn.M,pn.N);
            if dRow > 0
                pix(end-dRow+1:end,:) = 0;
                nb(1:dRow,:) = 0;
            elseif dRow < 0
                pix(1:-dRow,:) = 0;
                nb(end+dRow+1:end,:) = 0;
            end
            if dCol > 0
                pix(:,end-dCol+1:end) = 0;
                nb(:,1:dCol) = 0;
            elseif dCol < 0
                pix(:,1:-dCol) = 0;
                nb(:,end+dCol+1:end) = 0;
            end            
            plist = find(pix(:)');
            nlist = find(nb(:)');
            wlist = wt * ones(1,numel(plist));
            if ~isempty(pixMask)
                keep = pixMask(plist) & pixMask(nlist); %only keep edges where both plist and nlist are true
                plist = plist(keep);
                nlist = nlist(keep);
                wlist = wlist(keep);
            end
        end
        
        %This is similar to calcInds, except that rather than a single wt
        %per row/col offset, we can get an individual weight per element.
        %Find indices of neighbors
        %dRow, dCol: row and column offset defining a neighbor
        %Note: can also work if dRow and dCol are vectors of offsets
        %dRows, dCols: [1xQ] row and column offset of Q neighbors
        %wFull: [M x N] (if Q = 1) 
        %    OR {1 x Q} each element is [M x N], and wlist is found by
        %    using corresponding plist to look up values in this.  
        %plist: [1xP] indices of pixels having a neighbor
        %nlist: [1xP] indices of corresponding neighbors
        %wlist: [1xP] weights for each edge
        %pixMask: [] then all pixels in MxN
        %         [MxN] logical, then restricts pixels and neighbors to
        %         those both with pixMask(plist) == 1 and pixMask(nlist) == 1
        %Note: when using BK_SetNeighbors, dRow + M*dCol should be >=0 as
        %only the upper triangular  portion is considered
        function [plist, nlist, wlist] = calcIndsAll( pn, dRow, dCol, wFull, pixMask )
            if nargin<4
                if isscalar(dRow)
                    wFull = ones(pn.M, pn.N);
                else
                    wFull = cell(1,numel(dRow));
                end
            end;
            if nargin<5, pixMask = []; end;
            if numel(dRow)>1
                if numel(dRow) ~= numel(dCol) || numel(dRow) ~= numel(wFull)
                    error('dRow, dCol and wt must be vectors of equal length');
                end
                [pc,nc,wc] = cellfun( @(x,y,z) calcIndsAll(pn, x, y, z, pixMask), num2cell(dRow), num2cell(dCol), wFull, 'uniformoutput', false);
                plist = horzcat(pc{:});
                nlist = horzcat(nc{:});
                wlist = horzcat(wc{:});
                return;
            end
            %pix and nb will have the same size blocks of ones, just offset
            pix = ones(pn.M,pn.N);
            nb = ones(pn.M,pn.N);
            if dRow > 0
                pix(end-dRow+1:end,:) = 0;
                nb(1:dRow,:) = 0;
            elseif dRow < 0
                pix(1:-dRow,:) = 0;
                nb(end+dRow+1:end,:) = 0;
            end
            if dCol > 0
                pix(:,end-dCol+1:end) = 0;
                nb(:,1:dCol) = 0;
            elseif dCol < 0
                pix(:,1:-dCol) = 0;
                nb(:,end+dCol+1:end) = 0;
            end            
            plist = find(pix(:)');
            nlist = find(nb(:)');
            if iscell(wFull), wFull = wFull{1}; end;
            wlist = wFull(plist);
            if ~isempty(pixMask)
                keep = pixMask(plist) & pixMask(nlist); %only keep edges where both plist and nlist are true
                plist = plist(keep);
                nlist = nlist(keep);
                wlist = wlist(keep);
            end
        end
        
        %create a sparse array size [M*N x M*N] containing edges between
        %pixels
        %dRows, dCols: [1xQ] row and column offset of Q neighbors
        %wt: [1xQ] weight for each neighbor type, default: ones(1,Q)
        function sa = makeSparse( pn, dRows, dCols, wt, pixMask )
            if nargin<4 || isempty(wt)
                wt = ones(1,numel(dRows));
            end
            if nargin<5, pixMask = []; end;
            [plist,nlist,wlist] = pn.calcInds( dRows, dCols, wt, pixMask );
            sa = sparse(plist, nlist, wlist, pn.M*pn.N, pn.M*pn.N);
        end
        
        %rc: [rows ...; cols ...] non-zero elements of makeSparse
        function [rc, wts] = makeEdges( pn, varargin )
            [r,c,wts] = find( pn.makeSparse( varargin{:} ) );
            rc = [r(:)'; c(:)'];
            wts = wts(:)';
        end
        
        %Calculates pairwise potentials used by BK_SetPairwise
        %that are equivalent to using BK_SetNeighbors with makeSparse
        %pairwise: n x 6, each row: [i j e00 e01 e10 e11]
        %Here e01 and e10 are potentials (wt) for pixel i differing from
        %pixel j.
        function pairwise = getPairwise( pn, dRow, dCol, wt, pixMask )
            if nargin<4, wt = ones(1,numel(dRow)); end;
            if nargin<5, pixMask = []; end;
            [plist, nlist, wlist] = pn.calcInds( dRow, dCol, wt, pixMask );
            pairwise = [plist' nlist' zeros(numel(plist),1) wlist' wlist' zeros(numel(plist),1)];
        end
        
        function plot( pn, plist, nlist, sopt )
            if nargin < 4, sopt = {'r-'}; end;
            [py,px] = ind2sub([pn.M pn.N], plist);
            [ny,nx] = ind2sub([pn.M pn.N], nlist);
            ah = ishold;
            if ~ah, hold on; end
            plot([px;nx],[py;ny],sopt{:});
            if ~ah, hold off; end
        end
        
    end
    
end
