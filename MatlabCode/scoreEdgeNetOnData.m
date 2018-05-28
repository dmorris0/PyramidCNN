%scoreEdgeNetOnData
%
%  scoreEdgeNetOnData( imglist, seglistGT, edgefolder, edgeRadius, net, doGPU, redoEdges, redoSegs  )
%
%  Runs network on set of images. Outputs detected edges and runs
%  segmentation on them.  Uses GT to score them edges and segments.
%  Results are stored in text files in edgefolder.
%  This is a simplified version of applyEdgeNet.m, and is called by:
%  runEdgeSegNet.m
%  Use comparePerformance.m to display quantitative results
%
%  net:       trained dagNN network (MatConvNet)
%  imglist:   list of images
%  seglistGT:   list of gt segments
%  edgeRadius: radius around each GT edge pixel in which we seek best edge
%              score for precision/recall estimation
%  edgefolder: folder to save data
%  redoEdges, redoSegs: true (default) re-does CNN and/or segmentation even
%              if already done
%
%
%  See also: applyEdgeNet, runEdgeSegNet
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function scoreEdgeNetOnData( imglist, seglistGT, edgefolder, minfrac, edgeRadius, net, doGPU, redoEdges, redoSegs, crop, averageScales  )

if nargin < 4, minfrac = 0.1; end
if nargin < 5 || isempty(edgeRadius), edgeRadius = 2.5; end
if nargin < 6, net = []; end
if nargin < 7, doGPU = false; end
if nargin < 8, redoEdges = true; end  
if nargin < 9, redoSegs = true; end  
if nargin < 10, crop = []; end  %crop: [xmin xmax ymin ymax]
if ~isempty(crop) && ~numel(crop)==4, error('crop must be [] or 1x4'); end
if nargin < 11, averageScales = false; end

nonEdgeThresh = 0.98;

if ~iscell(imglist), error('Requires cell array of names'); end
nim = numel(imglist);

%% Use CNN to estimate leaf boundaries (or else load from file):
edgelist = runNetEdges( net, imglist, edgefolder, doGPU, redoEdges, crop, averageScales );

%% Do segmentation
if isempty(redoSegs)
    %then no segmentation, just edges
    seglist = [];
else
    seglist = runSegments( edgelist, minfrac, nonEdgeThresh, redoSegs );
end

%% Evaluation and plotting
f1 = figure(1); clf;
f2 = figure(2); clf;
f3 = figure(3); clf;
f4 = figure(4); clf;
f5 = figure(5); clf;
totpos = [];  %pos scores for edge pixels (max res edge)
totneg = [];  %neg scores for edge pixels
npixtp = 0;  %n pix true pos after segmentation
npixfp = 0;
npixfn = 0;
fprintf('\n==> Doing evaluation on %d edge and segment images\n', nim);
for its=1:nim
    fprintf('Reading %d of %d edge and segment images, ', its, nim);
    img = im2single( imread( imglist{its} ) );
    edge = im2single( imread( edgelist{its} ) );
    if size(edge,3)==3, edge = mean(edge,3); end
    if ~isempty(seglist)
        segLabels = colors2SegLabels( imread( seglist{its} ) );
    end
    
    [~,~,ext] = fileparts( seglistGT{its} );
    if strcmpi(ext,'.csv')
        segLabelsGT = csvread( seglistGT{its} );
    else
        segLabelsGT = imread( seglistGT{its} );
    end
    if ~isempty(crop)
        img = img(crop(3):crop(4),crop(1):crop(2),:);
        segLabelsGT = segLabelsGT(crop(3):crop(4),crop(1):crop(2),:);
    end

    figure(f1);
    clf
    im(img(:,:,1:3))
    if ~isempty(segLabelsGT), plotSegBoundaries( segLabelsGT ); end
    colorbar('off')
    title('Image with GT segments');
    
    figure(f2)
    clf
    [pos,neg] = segEdgePR( edge, segLabelsGT, edgeRadius);
    totpos = [totpos pos];  %keep track of all true postive and true negatives
    totneg = [totneg neg];
    APedge = prCurve( pos, neg );
    title(['Edge AP: ' num2str(APedge)]);

    if isempty(seglist)
        dice = 0;
        pixPrecision = 0;
        pixRecall = 0;
    else
        figure(f3)
        clf
        im(img(:,:,1:3))
        plotSegBoundaries(segLabels, {'m-','linewidth',2} );
        title('Segments on Image');
        
        figure(f4)
        [dice, pixPrecision, pixRecall, tp, fp, fn] = scoreEstSegs( segLabelsGT, segLabels, 1 - edge );
        npixtp = npixtp + sum(tp(:));
        npixfp = npixfp + sum(fp(:));
        npixfn = npixfn + sum(fn(:));
    end
    
    fprintf( '# Edge AP  \tedge rad \tminfrac \tDice     \tpixPr    \tpixRecall\n');
    fprintf( '%f \t%f \t%f \t%f \t%f \t%f\n', APedge, edgeRadius, minfrac, dice, pixPrecision, pixRecall);
        
    drawnow;
    
    fid = fopen(strrep(edgelist{its},'-edge.png','-AP.txt'),'w');
    fprintf(fid, '# Edge AP  \tedge rad \tminfrac \tDice     \tpixPr     \tpixRecall\n');
    fprintf(fid, '%f \t%f \t%f \t%f \t%f \t%f\n', APedge, edgeRadius, minfrac, dice, pixPrecision, pixRecall);
    fclose(fid);
end

allStatsFile = [edgefolder filesep 'AllStats.txt'];
%%
PR = npixtp / (npixtp + npixfp + 1e-5);
RL = npixtp / (npixtp + npixfn + 1e-5);
dice = 2*PR*RL / (PR + RL + 1e-5);
fid = fopen(allStatsFile,'w');
%figure(f5)
TAPedge = prCurve(totpos, totneg );
%title(['Cumulative Edge Average Precision: ' num2str(TAPedge,'%.3f')]);
fprintf(fid,'#TotEdge AP \tSegPixTP \tSegPixFP  \tSegPixFN  \tSegDice\tSegPrec\tSegRecall\n');
fprintf(fid,'%f  \t\t%d  \t%d  \t%d  \t%.3f  \t%.3f  \t%.3f\n', TAPedge, npixtp, npixfp, npixfn, dice, PR, RL );
fclose(fid);

%% Read output file and display result
%Illustrate how to do this:
fid = fopen(allStatsFile,'r');
vals = textscan(fid,'%f',7,'commentstyle','#');
fclose(fid);
vals = vals{1};
fprintf('#TotEdgePix AP\tSegPixTP \tSegPixFP  \tSegPixFN  \tSegDice\tSegPrec\tSegRecall\n');
fprintf('%f  \t\t%d  \t%d  \t%d  \t%.3f  \t%.3f  \t%.3f\n', vals );




end
