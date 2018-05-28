%runNetEdges
%
%  edgelist = runNetEdges( net, imglist, outfolder, doGPU, overwrite, crop, averageScales )
%
%  Run edge-detecting NN on list of images and saves predicted edges
%  net:  dagNN
%  imglist: cell array of image names
%  outfolder: folder to save edges to
%  averageScales: false (default) uses just highest resolution edge
%                 true: averages edges at all resolutions
%  edgeList:  names of edges
%
%  Daniel Morris, Oct 2017
%
%  See also: scoreEdgeNetOnData, runEdgeSegNet, runSegments
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function edgelist = runNetEdges( net, imglist, outfolder, doGPU, overwrite, crop, averageScales )

if nargin < 4, doGPU = false; end
if nargin < 5, overwrite = true; end
if nargin < 6, crop = []; end
if nargin < 7, averageScales = false; end

if ~iscell(imglist), error('Requires cell array of names'); end
nim = numel(imglist);
edgelist = cell(1,nim);
for its=1:nim
    base = getBaseName( outfolder, imglist{its}, its );
    edgelist{its} = [base '-edge.png'];
end
if isempty(net)
    fprintf('No network input, quitting\n');
    return
end
fprintf('\n==> Running trained network on %d images\n', nim);

%set predicted edges to store values:
[~,pinds] = setNetPredictions( net );

if doGPU
    net.move('gpu');
end

N = numel(pinds);
nupdown = max(findNDown( net ), N-1);

if ~isempty( outfolder )
    if exist( outfolder, 'dir' )
        fprintf('Will output to: %s\n', outfolder);
    else
        a = mkdir( outfolder );
        if ~a, error(['Cannot create: ' outfolder]); end
        fprintf('Creating folder for writing: %s\n', outfolder );
    end
end

for its=1:nim
    %skip files that have already been calculated
    if ~overwrite && exist(edgelist{its},'file')
        fprintf('Exists so skipping: %s\n', edgelist{its} );
        continue;
    end
    
    fprintf('Reading %d of %d images, ', its, nim);
    img = im2single( imread(imglist{its}) );

    if ~isempty(crop)
        warning('Cropping image');
        img = img(crop(3):crop(4),crop(1):crop(2),:);
    end
    %make sure image size is a multiple of a power of 2 needed to number of
    %half-sizes.  If not, then pad appropriately:
    osize = size(img);
    img = fixImDim( img, 2^nupdown );

    fprintf(' net eval on: %s', imglist{its});
    istart = tic;
    if doGPU
        net.eval( {'input', gpuArray( img )} );
        edge = pred2ProbMChan( gather( net.vars( pinds(1) ).value ) );
    else
        net.eval( {'input', img } );
        edge = pred2ProbMChan( net.vars( pinds(1) ).value );
    end
    if averageScales        
        for jts=2:N
            if doGPU
                nedge =  pred2ProbMChan( gather( net.vars( pinds(jts) ).value ) );
            else
                nedge =  pred2ProbMChan( net.vars( pinds(jts) ).value );
            end
            edge = edge + imresize( nedge, 2^(jts-1), 'bilinear' );
        end
        edge = edge/N;
    end
    
    fprintf(', done in: %5.1f sec. ', toc(istart) );
    
    %revert to size of image:
    edge = edge(1:osize(1), 1:osize(2), : );
    
    imwrite( uint16(edge*(2^16-1)),edgelist{its});
    fprintf('Saved edges\n');
end
    
if doGPU
    net.move('cpu');
end


end

function base = getBaseName( outfolder, imname, its )
if ischar( imname )
    [~,basen] = fileparts( imname );
else
    basen = ['img-' num2str(its,'%03d')];
end
base = [outfolder filesep basen];
end

    
function n = findNDown( net )

%find number of times network downsamples an image:
%This assumes a downsample layer name contains: '_dn'
n = 0;
for its=1:numel(net.layers)
    if contains( net.layers(its).name, '_dn')
        n = n+1;
    end
end

end
