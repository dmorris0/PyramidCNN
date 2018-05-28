%prCurve
%
%  [AP,h] = prCurve(posData, negData, sval)
%  [AP,h] = prCurve(posData, negData, sval, nPos)
%
%  Plots precision versus recall
%  posData, negData: scores for positives and negatives
%  sval:      cell array for plot, default: {'b-'}
%             if {} then does not plot and only returns AP
%  nPos:      if posData do not cover all ground truth positives, then
%             specify manually the total number of positives, default:
%             numel(posdata)
%  precision: TP/(TP + FP)
%  recall:    TP/(TP + FN) = TP / nPos
%  AP:        average precision
%  h:         handle for plot
%
%  Daniel Morris, Aug 2013
%
%  See also: rocCurve
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function [AP, h] = prCurve(posData, negData, sval, nPos)

if nargin<4, nPos = numel(posData); end
if nargin<3, sval = {'b-'}; end
if ischar(sval), sval = {sval}; end

%eliminate pos data with -inf (after we calculate nPos above)
posData = posData(posData>-inf);

data = [posData(:)' negData(:)'];
[scores,order] = sort(data,'descend');
TP = [ones(1,numel(posData)) zeros(1,numel(negData))];
TP = TP(order);
FP = 1-TP;
%Here's an issue.  If scores do not change between successive elements, it
%is wrong to treat FP and TP counts incrementally.  This produces
%non-monotonically decreasing PR curves. Rather we need to group
%all counts with the same score together:
deltas = (scores - scores([2:end end]) ) > 0;
cumdeltas = cumsum(deltas) + 1;
nd = max(cumdeltas);
np = numel(data);
G = sparse(cumdeltas, 1:np, ones(1,np), nd, np );
TP = TP * G';
FP = FP * G';

cumTP = cumsum(TP);
cumFP = cumsum(FP);
detPos = max(cumTP + cumFP, 1);  %avoid divide by 0

precision = cumTP ./ detPos;
recall = cumTP / nPos;

%The curve can be extended out horizontally to the left:
if recall(1) > 0
    precision = precision([1 1:end]);
    recall = [0 recall];
end

%calculate average precision by summing area under curve:
doLowerBound = false;
if doLowerBound
    barheight = precision(2:end);
else
    barheight = 0.5*(precision(1:end-1)+precision(2:end));
end
barwidth = recall(2:end) - recall(1:end-1);
AP = barheight * barwidth';

if ~isempty(sval)
    h = plot(recall, precision, sval{:});
    xlabel('Recall');
    ylabel('Precision');
end

end



