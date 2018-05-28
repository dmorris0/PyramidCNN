%runLeafNet
%
%  Script to illustrate running a pre-trained network on data
%  Also uses boundary output to do instance segmentation, and evaluates
%  result using ground truth segments
%
%
%  Daniel Morris, May 2017
%
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License

%%
topFolder = '..';  %assumes we are running from Scripts folder
dataPath = [topFolder filesep 'Data'];
dataName = 'DenseLeaves';
testFolder = [dataPath filesep dataName filesep 'test'];
modelPath = [topFolder filesep 'TrainedModels'];
outputPath = [topFolder filesep 'Output'];
modelName = 'net_90_183';

%% Check that dataset is downloaded and if not download it
if ~exist([dataPath filesep dataName],'dir')
    fprintf('Dataset not found, downloading %s ... ', dataName);
    cpwd = pwd;
    cd(dataPath)
    unzip(['https://www.egr.msu.edu/denseleaves/Data/' dataName '.zip']);
    cd(cpwd);
    fprintf('done!\n');
end

%% Specify data to run net on:
plantNames = getImageList( testFolder, '_img','png' );
segNames = cellfun(@(x)strrep(x,'_img','_seg'),plantNames,'uniformoutput',false);

fprintf('Will run on %d images in folder: %s\n', numel(plantNames), testFolder );

%% Load network
netName = [modelPath filesep modelName];
a = load(netName, 'net') ;
net = dagnn.DagNN.loadobj(a.net) ;
fprintf('Loaded net: %s\n', netName );

%% Run network on images, do segmentation and evaluate results:
if ~exist(outputPath,'dir'), mkdir(outputPath); end
outfolder = [outputPath filesep modelName];
minfrac = 0.1;    %Used in merging segments to build leaf estimates
edgeRadius = 2.5; %Used in evaluation of edge boundary pixels
doGPU = false;    %if have GPU and compiled MatConvNet with GPU set to true
redoEdges = false;
redoSegs = false;
doAverage = false;
scoreEdgeNetOnData( plantNames, segNames, outfolder, minfrac, edgeRadius, net, doGPU, redoEdges, redoSegs , [], doAverage  );

