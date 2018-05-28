%mergeSegs
%
%  nsegim = mergeSegs( segim, density, minfrac, method=='larger' )
%  nsegim = mergeSegs( segim, distim, maxdist, method=='maxdist' )
%
%  Merge segments based on distim
%  segim:  [M x N] segment labels
%  distim: [M x N] distance from edge for each pixel
%  maxdist: 1x1 if distance from edge along a boundary is greater than this
%               then will merge segments
%  method:   'maxdist': does max dist method
%          'larger':  does larger seg merge method
%
%  Daniel Morris, Jun 2017
%
%  See also: showNetWeights, showNetActivation, showNetPRcurve, showNetSegs
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function segim = mergeSegs( segim, density, minfrac, method, fignum )

if nargin < 4, method = 'larger'; end
if nargin < 5, fignum = []; end

Ns = max(segim(:));
if Ns==numel(segim)
    warning('One segment per pixel -- something is wrong');
end
[M,N] = size(segim);

pn = pixNeighbor(M, N);
[plist,qlist] = pn.calcInds( [-1 0 1 1], [1 1 1 0] );  %8-connected

switch method
    
    case 'maxdist'
        
        segim = maxConMerge( segim, density, minfrac, Ns, plist, qlist, fignum );
        
    case 'larger'
        
        segim = largerMerge( segim, density, minfrac, Ns, plist, qlist, 10, fignum );
        
    otherwise
        error('Invalid mopt');
end

end

function fignum = plotDenSeg( density, segim, lowind, highind, fignum)
if ~isempty(fignum)
    figure(fignum)
    clf
    imc = imCmap(density, 2);
    smerge = ismember( segim, [lowind; highind] );
    imc.setPixBand( smerge, 2);  %segments to merge
    imc.show;
    plotSegBoundaries( segim );
    fignum = fignum+1;
    [u,v] = meshgrid(1:size(segim,2), 1:size(segim,1));
    hold on
    for its=1:numel(lowind)
        xl = mean( u( segim==lowind(its) ) );
        yl = mean( v( segim==lowind(its) ) );
        xh = mean( u( segim==highind(its) ) );
        yh = mean( v( segim==highind(its) ) );
        
        plot([xl;xh],[yl;yh],'m-');
    end
    hold off
    drawnow
end
end

function [segim, nseglist, fignum] = largerMerge( segim, edensity, minfrac, Ns, plist, qlist, nrepeat, fignum )

nseglist(1) = Ns;
fprintf('Merging segs with minfrac: %.2f, Ns: %d', minfrac, Ns);
for its=1:nrepeat
    %find which edges connect segments:
    connected = segim(plist)~=segim(qlist);
    
    %store total connection between pairs of segments
    %where a connection is the minimum of edensity at either end of the edge
    sumCon = accumarray([segim(plist(connected)') segim(qlist(connected)')], ...
        min(edensity(plist(connected)'), edensity(qlist(connected)') ), ...
        [Ns Ns], ...
        @sum);
    
    %get both-way connection betweeh segments:
    segCon = sumCon + sumCon';
    
    %now find total number of connections from each segment
    totCon = accumarray([segim(plist(connected)');segim(qlist(connected)')], ...
        ones(2*sum(connected(:)),1), ...
        [Ns 1] );
    
    %fraction of edges connecting to neighbors:
    fracCon = bsxfun(@rdivide, segCon, totCon );
    
    %find number of pixels per segment:
    npix = accumarray(segim(:), ones(numel(segim(:)),1), [Ns,1] );
    [~,sind] = sort(npix);
    
    %sorted fracCon and eliminate connections to neighbors with fewer pixels
    sfracCon = triu( fracCon(sind,sind) );
    
    %Get maximum fraction neighbor (with greater number of pixels):
    [mfrac,mind] = max( sfracCon, [], 2);
    
    %do plotting
    fignum = plotDenSeg( edensity, segim, sind(mfrac > minfrac), sind(mind(mfrac > minfrac)), fignum );    
    
    %now do connected components
    G = sparse(sind(mfrac > minfrac), sind(mind(mfrac > minfrac)), 1, Ns, Ns);
    [~,segLabels] = graphconncomp(G+G','Directed',false);  %make symmstric
    
    %Now give pixels from merged segments the same numbers:
    nsegim = zeros(size(segim));
    L = unique(segLabels);
    %New number of segments:  
    Ns = numel(L);
    fprintf(', %d',Ns);
    nseglist(end+1) = Ns;
    for jts=1:Ns
        nsegim( ismember(segim, find(segLabels==jts) ) ) = jts;
    end
    segim = nsegim;
    if nseglist(end)==nseglist(end-1)
        %done;
        break;
    end
end
fprintf('\n');

end



function segim = maxConMerge( segim, distim, maxdist, Ns, plist, qlist, fignum )

%find which edges connect segments:
connected = segim(plist)~=segim(qlist);

%Note: Sometimes I get a memory error here becaues Ns is equal to the
%number of pixels in segim.  That is too many segments!
%store maximum connection between segments
%where a connection is the minimum of distim at either end of the edge
maxCon = accumarray([segim(plist(connected)') segim(qlist(connected)')], ...
    min(distim(plist(connected)'), distim(qlist(connected)') ), ...
    [Ns Ns], ...
    @max);

%make symmetric and keep strongest edge connection between each segment:
maxCon = max(maxCon,maxCon');

%now merge rows/columns
ind = 1;
while ind < size(maxCon,1)
    %consider segment ind (starting at lowest and working up)
    %find which segments this connects to strongly enough:
    ctomerge = maxCon(ind,:) > maxdist;
    if sum(ctomerge) > 0
        %merge all segments with strong connections
        %rename their pixels with segs(ind):
        segim( ismember(segim, segs(ctomerge)) ) = segs(ind);
        %add all their connections to ind:
        maxCon(ind,:) = max( [maxCon(ind,:);maxCon(ctomerge,:)], [], 1);
        maxCon(:,ind) = maxCon(ind,:)';
        %delete their connections
        maxCon(ctomerge,:) = 0;
        maxCon(:,ctomerge) = 0;
        %make sure ind does not connect to itself:
        maxCon(ind,ind) = 0;
    else
        %if ind has no more connections, go onto next
        ind = ind+1;
    end
end

end

