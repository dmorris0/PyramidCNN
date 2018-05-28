%scoreEstSegs
%
%  [dice, precision, recall, tp, fp, fn] = scoreEstSegs( gtSegs, estSegs )
%  [dice, precision, recall, tp, fp, fn] = scoreEstSegs( gtSegs, estSegs, imForPlot )
%
%  Scores a proposed segmentation image.  Finds association between ground
%  truth segments and estimated segments that maximizes the total segment
%  dice scores.  Given this association each pixel is labeled as a tp, fp
%  or fn.  Using this precision, recall and dice are calculated for the
%  image.
%  gtSegs: [M x N] segment labels, 1:Ns segments (zero values are ignored)
%  estSegs: [M x N] estimated segment labels, values between 1 and Ne
%  imForPlot:  [MxN] optional: if input then plots this image with the
%                   gtSegs and corresponding estSegs 
%  diceVal: 2*Precision*Recall/(Precision+Recall)
%  tp, fp, fn: [M x N] indicate if pixels are true positives, false
%              positives of false negatives, given their association with
%              gtSegs
%
%
%  Daniel Morris, Jan 2018
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function [dice, precision, recall, tp, fp, fn] = scoreEstSegs( gtSegs, estSegs, imForPlot, doQuarter )

if nargin < 3, imForPlot = []; end
if nargin < 4, doQuarter = false; end

%sometimes estSegs will be larger (as images dimension need to be powers of
%2).  So crop it to gtSegs size:
if numel(gtSegs) < numel(estSegs)
    estSegs = estSegs( 1:size(gtSegs,1), 1:size(gtSegs,2) );
    if ~isempty(imForPlot)
        imForPlot = imForPlot( 1:size(gtSegs,1), 1:size(gtSegs,2), : );
    end
end

%Create association between gtSegs and estSegs:
[diceVals, segAssoc] = maxSegDice( gtSegs, estSegs );

tp = false(size(gtSegs));
fp = false(size(gtSegs));
fn = false(size(gtSegs));
for its=1:numel(diceVals)
    if diceVals(its) > 0
        tp = tp | (gtSegs == segAssoc(its,1) & estSegs == segAssoc(its,2));
        fp = fp | (gtSegs ~= segAssoc(its,1) & estSegs == segAssoc(its,2));
        fn = fn | (gtSegs == segAssoc(its,1) & estSegs ~= segAssoc(its,2));
    else
        %no positives, only negatives:
        fn = fn | (gtSegs == segAssoc(its,1));
    end
end

precision = sum(tp(:)) / (sum(tp(:)) + sum(fp(:)));
recall = sum(tp(:)) / (sum(tp(:)) + sum(fn(:)));
dice = 2*precision*recall / (precision + recall);

if ~isempty(imForPlot)
    cmap = [0 1 0; %tp: green
        1 0 0;     %fn (misses): red
        0 0 1];    %fp (mistakes): blue
    
    if doQuarter
        sim = smoothim( imForPlot, 2);
        imc = imCmap( sim(2:4:end,2:4:end), 4, 0.4, cmap );
        imc.setPixBand( tp(2:4:end,2:4:end), 2 );
        imc.setPixBand( fn(2:4:end,2:4:end), 3 );
        imc.setPixBand( fp(2:4:end,2:4:end), 4 );
        imc.show;
        axis off
        
        keep = ismember( estSegs, segAssoc(:,2) );
        plotSegBoundaries( estSegs(2:4:end,2:4:end).*keep(2:4:end,2:4:end), {'m-'}, false );
        title(['Quarter-size, Dice: ' num2str( dice, ' %.3f' ) ', P: ' num2str(precision,'%.2f') ', R: ' num2str(recall,'%.2f') ] );
    else
        imc = imCmap( imForPlot, 4, 0.4, cmap );
        imc.setPixBand( tp, 2 );
        imc.setPixBand( fn, 3 );
        imc.setPixBand( fp, 4 );
        imc.show;
        axis off
        
        keep = ismember( estSegs, segAssoc(:,2) );
        plotSegBoundaries( estSegs.*keep, {'m-'}, false );
        title(['Dice: ' num2str( dice, ' %.3f' ) ', P: ' num2str(precision,'%.2f') ', R: ' num2str(recall,'%.2f') ] );
    end
end


end



