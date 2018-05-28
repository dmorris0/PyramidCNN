%setNetPredictions
%
%  [activations, pinds] = setNetPredictions( net, setPrecious, predname )
%  [activations, pinds] = setNetPredictions( net, setPrecious, activations )
%
%  Sets prediction activations to precious = 1, also returns activation
%  names and indices.  This way running the network will save the
%  predictions, and they can be accessed with: net.vars(pinds(its)).value
%  setPrecious: true (default) sets predictions to true 
%  predname:    base name of prediction 'predlabel' is default.  
%               identifies prediction with: strfind( net.vars(its).name, predname )
%  activations: cell array of names of predictions
%  pinds:       indices into net.vars for each prediction
%
%  Note: to set network to remember all activations use:     
%        net.conserveMemory = false;
%
%  Daniel Morris, Nov 2017
%
%  See also: getNetVars, plotNetConnections, runNetFolder, removeFinalLayers
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function [activations, pinds] = setNetPredictions( net, setPrecious, predname )

if nargin < 2, setPrecious = true; end
if nargin < 3, predname = 'predlabel'; end

activations = {};
pinds = [];
if iscell(predname)
    activations = predname;
    for its=1:numel(predname)
        pinds(its) = net.getVarIndex(activations{its});
        net.vars( pinds(its) ).precious = 1;
    end
else
    pinds = getAllPred( net, predname );
    activations = cell(1,numel(pinds));
    for its=1:numel(pinds)
        activations{its} = net.getVar(pinds(its)).name;
        net.vars( pinds(its) ).precious = 1;
    end
end
%     inc = 0;
%     while true
%         inc = inc+1;
%         ind = net.getVarIndex([predname num2str(inc)]);
%         if isnan(ind)
%             break;
%         else
%             pinds(end+1) = ind;
%             activations{end+1} = net.getVar(ind).name;
%         end
%     end

if setPrecious
    fprintf('=========================\n');
    fprintf('Set %d layers precious=1: ',numel(pinds));
    fprintf('%s ',activations{:});
    fprintf('\n');
end

end


function ind = getAllPred( net, basename )

ind = [];
for its=1:numel(net.vars)
    if ~isempty(strfind( net.vars(its).name, basename ))
        ind(end+1) = its;
    end
end

end

