%imageTrainData
%
%  Class for storing and accessing image and label info for training a CNN
%
%  itd = imageTrainData
%
%  itd.init( imdb )
%
%  
%
%  Daniel Morris, Nov 2017
%
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
classdef imageTrainData < handle
    
    properties
        id               %1xQ data id
        data             %[MxNxPxQ] image data
        data2            %[MxNxRxQ] additional data besides images
        set              %1xQ 1 = train, 2 = val
        instanceWeights  %[MxNx1xQ] pixel weights
        labelData        %{nlevels x nPerLevel} cell array of labels
        labelNames       %{nlevels x nPerLevel} Names of labels
        labelWeights     %{nlevels x nPerLevel} weights for each label at each level.  These can either be per category or per pixel, see itd.calcWeights for per category        
        labelWeightNames %{nlevels x nPerLevel} Names of label weights (when weights stored as a per-label-pixel weight)
        reweights        %{nlevels x nperLevel} scale factor to multiply weights by to give more weight to hard negatives
        opts             %cell array of options used in defining data
    end
    
    methods
        %constructor
        function itd = imageTrainData( varargin )
            if nargin > 0
                itd.init( varargin{:} );
            end
        end
        
        function init(itd, imdb )
            if ischar( imdb )
                itd.load( imdb );
            elseif isstruct(imdb)
                itd.id = imdb.images.id;
                itd.data = imdb.images.data;
                itd.data2 = [];
                itd.set = imdb.images.set;
                itd.instanceWeights = imdb.images.instanceWeights;
                itd.labelData = {imdb.images.label};
                itd.labelNames = {'label'};
                itd.labelWeights = {};
                itd.labelWeightNames = {};
                if isfield( imdb,'myOpts')
                    itd.opts = imdb.myOpts;
                end
            else
                error('Unknown data type input');
            end
        end
        
        function reset(itd)
            itd.id = [];
            itd.data = [];
            itd.data2 = [];
            itd.set = [];
            itd.instanceWeights = [];
            itd.labelData = {};
            itd.labelNames = {};
            itd.labelWeights = {};
            itd.opts = struct;
        end
        
        function nl = nlevels(itd)
            nl = size(itd.labelData,1);
        end
        
        function nperl = nperlevel(itd)
            nperl = size(itd.labelData,2);
        end
        
        %for backward compatibility with imdb:
        function firstlabel = label( itd )
            firstlabel = itd.labelData{1};
        end
        
        function n = ntrain( itd )
            n = sum(itd.set==1);
        end
        
        function n = nval( itd )
            n = sum(itd.set==2);
        end
        
        %trainOrVal: 'train' (default) or 'val'
        function img = getIm( itd, ind, trainOrVal, channels )
            if strcmpi(trainOrVal,'train')
                ilist = find(itd.set==1,ind);
            else
                ilist = find(itd.set==2,ind);
            end
            if numel(ilist) < ind
                img = [];
            else
                if nargin < 4
                    img = itd.data(:,:,:,ilist(end));
                else
                    img = itd.data(:,:,channels,ilist(end));
                end
            end
        end
        
        %find n'th element of train or val set
        function ind = getInd( itd, ind, trainOrVal )
            if nargin < 3, trainOrVal = 'train'; end
            if strcmpi(trainOrVal,'train')
                ilist = find(itd.set==1);
            else
                ilist = find(itd.set==2);
            end
            ind = ilist(ind);
        end
        
        function save(itd, filename)
            save(filename,'itd','-v7.3');
        end
        
        function load(itd, filename)
            a = load(filename);
            if ~isfield(a,'itd')
                error(['Did not find imageTrainData in: ', filename]);
            end
            itd.copy(a.itd);
        end
        
        function copy(itd, oitd)
            fn = properties(itd);
            for its = 1:numel(fn)
                itd.(fn{its}) = oitd.(fn{its});
            end
        end
        
        function mem = getMem( itd )
            plist = properties( itd );
            mem = 0;
            for its=1:numel(plist)
                curprop = itd.(plist{its});
                s = whos('curprop');
                mem = mem + s.bytes;
            end
            if nargout==0
                if mem > 1e6
                    fprintf('Data memory use: %.2f Gb\n', mem/1e9 );
                else
                    fprintf('Data memory use: %.2f Mb\n', mem/1e6 );
                end
                clear mem
            end
        end                
        
        %This stores a category weight in labelWeights.  It is also
        %possible to store a per-pixel weight
        function calcWeights( itd )
            for its = 1:size(itd.labelData,1)
                for jts = 1:size(itd.labelData,2)
                    itd.labelWeights{its,jts} = calcClassWeights( itd.labelData{its,jts}(:,:,:,itd.set==1) );
                end
            end
        end
        
        function data = applyDropout( ~, data, mask, dropchan1, dropchan2, dropChan2Sign )
            if nargin < 6, dropChan2Sign = 1; end
            if ~isempty(mask) 
                if (dropchan1~=0) 
                    %set mask pixels in channel dropchan0 to 0
                    layer = data(:,:,dropchan1,:);
                    layer(mask) = 0;
                    data(:,:,dropchan1,:) = layer;
                end
                if (dropchan2~=0) 
                    %set mask pixels in channel dropchan1 to 1
                    layer = data(:,:,dropchan2,:);
                    layer(mask) = 1;
                    data(:,:,dropchan2,:) = dropChan2Sign * layer;
                    %if dropChan1Sign is -1, this makes second channel span
                    %from 0 to -1, so max values are closer points
                end
            end
        end
        
        %given losses over data calculate reweighting to up-weight hard
        %negatives
        %labels: [1 2] to reweight
        %hardfrac: fraction of hard negatives, default 0.1
        %hardweight: default 5 (rewight them by this).  
        % Note hardweight * hardfrac must be < 1 
        function setReweights( itd, epochloss, labels )
            if nargin < 3, labels = [1 2]; end
            if itd.opts.hardNegFrac * itd.opts.hardNegWeight >= 1
                error('hardNegFrac * hardNegWeightmust be < 1');
            end
            %weight on easy points (so that total weighted sum is same)
            %easyweight*(1-hardfrac) + hardweight*hardfrac = 1
            easyweight = (1-itd.opts.hardNegWeight*itd.opts.hardNegFrac) / (1-itd.opts.hardNegFrac);
            for its=1:numel(itd.labelWeights)
                if isa(epochloss{1},'gpuArray')
                    itd.reweights{its} = gpuArray( ones(size(itd.labelData{its}),'single') );
                else
                    itd.reweights{its} = ones(size(itd.labelData{its}),'single');
                end
                if ischar(labels)  %this is for regression, all pixels treated the same
                    pix = itd.labelWeights{its} > 0;
                    %find quantile for hardfrac:
                    thresh = quantile(epochloss{its}(pix), 1-itd.opts.hardNegFrac);
                    %set reweights as either easyweight or hardweight
                    %depending on their loss:
                    itd.reweights{its}(pix) = single(epochloss{its}(pix) < thresh) * easyweight + ...
                        single(epochloss{its}(pix) >= thresh) * itd.opts.hardNegWeight;
                else
                    for jts=1:numel(labels)
                        %find pixels for a given label:
                        pix = itd.labelData{its}==labels(jts) & itd.labelWeights{its} > 0;
                        %find quantile for hardfrac:
                        thresh = quantile(epochloss{its}(pix), 1-itd.opts.hardNegFrac);
                        %set reweights as either easyweight or hardweight
                        %depending on their loss:
                        itd.reweights{its}(pix) = single(epochloss{its}(pix) < thresh) * easyweight + ...
                            single(epochloss{its}(pix) >= thresh) * itd.opts.hardNegWeight;
                    end
                end
            end
            fprintf('Set itd.reweights using epochloss\n');
        end
        
        function [goodweights, iter] = loadLatestReweights(itd, folder )

            [iter,second] = findLastCheckpoint( folder );
            try
                goodweights = itd.loadReweights( folder, iter );
            catch
                %sometimes there is an error in saving reweights (due to
                %running out of time and the process being killed).  In
                %thta case load the previous reweights
                warning(['Unable to load latest reweights from iter: ' num2str(iter) ', trying from iter: ' num2str(second)]);
                iter = second;
                goodweights = itd.loadReweights( folder, iter );                
            end            
            
            function [epoch,secondToLast] = findLastCheckpoint(folder)
                list = dir(fullfile(folder, 'reweights-*.mat')) ;
                tokens = regexp({list.name}, 'reweights-([\d]+).mat', 'tokens') ;
                epoch = cellfun(@(x) sscanf(x{1}{1}, '%d'), tokens) ;
                epoch = sort([epoch 0],'descend');
                secondToLast = epoch(min(2,numel(epoch)));
                epoch = epoch(1);
            end

        end
        
        function goodweights = loadReweights(itd, folder, iter)
            filename = [folder filesep 'reweights-' num2str(iter,'%02d') '.mat'];
            if exist(filename,'file')
                a = load(filename);
                goodweights = true;
                if numel(a.reweights) ~= numel(itd.labelWeights)
                    goodweights = false;
                else
                    for its=1:numel(a.reweights)
                        if numel(a.reweights{its}) ~= numel(itd.labelWeights{its})
                            goodweights = false;
                            break;
                        end
                    end
                end
                if goodweights
                    if isa(itd.labelWeights{1},'gpuArray')
                        itd.reweights = moveToGPU( a.reweights );
                    else
                        itd.reweights = a.reweights;
                    end
                    fprintf('Loaded and initialized reweights from: %s\n', filename );
                else
                    warning(['Reweights invalid: ' filename]);
                end
            else
                goodweights = false;
            end
        end
        
        function saveReweights(itd, folder, iter)
            if ~isempty(itd.reweights) 
                if isa(itd.reweights{1},'gpuArray')
                    reweights = gather( itd.reweights );
                else
                    reweights = itd.reweights;
                end
                filename = [folder filesep 'reweights-' num2str(iter,'%02d') '.mat'];
                save(filename,'reweights','-v7.3');
                fprintf('Saved: %s\n', filename);
            end
        end

        %return required input to net.eval() for particular batch of data
        %can provide.  See plotNetIO.m
        function binput = itdBatch(itd, batch, gpus )
            if nargin < 3, gpus = []; end
            images = itd.data(:,:,:,batch);
            if itd.opts.dropFrac > 0 && ...
                isfield(itd.opts, 'dropChanList') && sum(itd.opts.dropChanList) > 0
                %mask = 1 for dropped pixels
                mask = rand(size(images,1), size(images,2), 1, size(images,4) ) < itd.opts.dropFrac;
            else
                mask = [];
            end
            images = itd.applyDropout( images, mask, itd.opts.dropChanList(1), itd.opts.dropChanList(2), itd.opts.dropChan2Sign );            
            if ~isempty(gpus)
                images = gpuArray( images );
            end           
            %first allocate binput:
            nlabels = numel(itd.labelNames);
            nlabelweights = numel(itd.labelWeightNames);
            nadd = 0;
            if ~isempty(itd.instanceWeights)
                nadd = nadd + 1;
            end
            if ~isempty(itd.data2)
                nadd = nadd + 1;
            end
            binput = cell(1, 2*(nlabels + nadd + nlabelweights) );
            %now fill binput:
            binput(1:2) = {'input', images};            
            inc = 2;
            if ~isempty(itd.instanceWeights)
                if isempty(gpus)                    
                    binput(inc+[1 2]) = {'instanceWeights', itd.instanceWeights(:,:,:,batch)};
                else
                    binput(inc+[1 2]) = {'instanceWeights', gpuArray( itd.instanceWeights(:,:,:,batch) )};
                end
                inc = inc+2;
            end
            if ~isempty(itd.data2)
                images2 = itd.data2(:,:,:,batch);
                images2 = itd.applyDropout( images2, mask, itd.opts.dropChanList(3), itd.opts.dropChanList(4), itd.opts.dropChan2Sign  );
                if nargin >= 3 && ~isempty(gpus)
                    images2 = gpuArray( images2 );
                end
                binput(inc+[1 2]) = {'data2', images2};
                inc = inc+2;
            end
            for its=1:nlabels
                if isempty(gpus)
                    binput(inc+[1 2]) = {itd.labelNames{its}, itd.labelData{its}(:,:,:,batch)};
                else
                    binput(inc+[1 2]) = {itd.labelNames{its}, gpuArray( itd.labelData{its}(:,:,:,batch)) };
                end
                inc = inc+2;
            end
            for its=1:nlabelweights
                %If weights have a single value per pixel but multiple
                %labels, then expand them to have 1 value per label
                n3 = size(itd.labelData{its},3);
                if isempty(gpus)              
                    if isempty(itd.reweights)
                        if size(itd.labelWeights{its},3) == 1 && n3 >1
                            binput(inc+[1 2]) = {itd.labelWeightNames{its}, itd.labelWeights{its}(:,:,ones(1,n3),batch)};
                        else
                            binput(inc+[1 2]) = {itd.labelWeightNames{its}, itd.labelWeights{its}(:,:,:,batch)};
                        end
                    else
                        binput(inc+[1 2]) = {itd.labelWeightNames{its}, itd.labelWeights{its}(:,:,:,batch).*itd.reweights{its}(:,:,:,batch) };
                    end
                else
                    if isempty(itd.reweights)
                        if size(itd.labelWeights{its},3) == 1 && n3 >1
                            binput(inc+[1 2]) = {itd.labelWeightNames{its}, gpuArray( itd.labelWeights{its}(:,:,ones(1,n3),batch)) };
                        else
                            binput(inc+[1 2]) = {itd.labelWeightNames{its}, gpuArray( itd.labelWeights{its}(:,:,:,batch)) };
                        end
                    else
                        binput(inc+[1 2]) = {itd.labelWeightNames{its}, gpuArray( itd.labelWeights{its}(:,:,:,batch).*itd.reweights{its}(:,:,:,batch) ) };
                    end
                end
                inc = inc+2;
            end
        end   
        
        %return required input to net.eval() for particular batch of data
        %can provide.  See plotNetIO.m
        %Got rid of gpus option since now use moveToGPU for itd -- actually
        %this only works when the data is not too large
        function binput = itdBatchGPU(itd, batch )
            images = itd.data(:,:,:,batch);
            if itd.opts.dropFrac > 0 && ...
                isfield(itd.opts, 'dropChanList') && sum(itd.opts.dropChanList) > 0
                %mask = 1 for dropped pixels
                mask = rand(size(images,1), size(images,2), 1, size(images,4) ) < itd.opts.dropFrac;
            else
                mask = [];
            end
            images = itd.applyDropout( images, mask, itd.opts.dropChanList(1), itd.opts.dropChanList(2) );            
            %first allocate binput:
            nlabels = numel(itd.labelNames);
            nlabelweights = numel(itd.labelWeightNames);
            nadd = 0;
            if ~isempty(itd.instanceWeights)
                nadd = nadd + 1;
            end
            if ~isempty(itd.data2)
                nadd = nadd + 1;
            end
            binput = cell(1, 2*(nlabels + nadd + nlabelweights) );
            %now fill binput:
            binput(1:2) = {'input', images};            
            inc = 2;
            if ~isempty(itd.instanceWeights)
                binput(inc+[1 2]) = {'instanceWeights', itd.instanceWeights(:,:,:,batch)};
                inc = inc+2;
            end
            if ~isempty(itd.data2)
                images2 = itd.data2(:,:,:,batch);
                images2 = itd.applyDropout( images2, mask, itd.opts.dropChanList(3), itd.opts.dropChanList(4) );
                binput(inc+[1 2]) = {'data2', images2};
                inc = inc+2;
            end
            for its=1:nlabels
                binput(inc+[1 2]) = {itd.labelNames{its}, itd.labelData{its}(:,:,:,batch)};
                inc = inc+2;
            end
            for its=1:nlabelweights
                if isempty(itd.reweights)
                    binput(inc+[1 2]) = {itd.labelWeightNames{its}, itd.labelWeights{its}(:,:,:,batch)};
                else
                    binput(inc+[1 2]) = {itd.labelWeightNames{its}, itd.labelWeights{its}(:,:,:,batch).*itd.reweights{its}(:,:,:,batch) };
                end
                inc = inc+2;
            end
        end            

        %create very simple image data
        function initDemo( itd, sel )
            if nargin < 2, sel = 1; end
            itd.reset;
            switch sel
                case 1
                    %1 level, 1 label
                    itd.data = cat(4, eye(8,'single'), triu( ones(8,'single') ) );
                    itd.set = [1 2];
                    itd.instanceWeights = [];
                    itd.labelData = cell(1);
                    itd.labelData{1} = cat(4, single(2)-eye(8,'single'), single(2)-triu( ones(8,'single') ) );
                    itd.labelNames = {'L1'};
                case 2
                    %1 level 2 labels
                    itd.data = cat(4, eye(8,'single'), triu( ones(8,'single') ) );
                    itd.set = [1 2];
                    itd.instanceWeights = [];
                    itd.labelData = cell(1,2);
                    itd.labelData{1} = cat(4, single(2)-eye(8,'single'), single(2)-triu( ones(8,'single') ) );
                    itd.labelData{2} = cat(4, single(2)-eye(8,'single'), single(2)-triu( ones(8,'single') ) );
                    itd.labelNames = {'L1','L2'};
                case 3
                    %2 levels, each 1 label
                    itd.data = cat(4, eye(8,'single'), triu( ones(8,'single') ) );
                    itd.set = [1 2];
                    itd.instanceWeights = [];
                    itd.labelData = cell(2,1);
                    itd.labelData{1} = cat(4, single(2)-eye(8,'single'), single(2)-triu( ones(8,'single') ) );
                    itd.labelData{2} = cat(4, single(2)-eye(8,'single'), single(2)-triu( ones(8,'single') ) );
                    itd.labelNames = {'L1';'L2'};
                case 4
                    itd.id = 1:2;
                    itd.data = repmat( cat(4, eye(8,'single'), triu( ones(8,'single') ) ), [1 1 1 5]);
                    itd.set = [ones(1,8) 2*ones(1,2)];
                    itd.instanceWeights = [];
                    itd.labelData = cell(2,1);
                    itd.labelData{1} = repmat( cat(4, single(2)-eye(8,'single'), single(2)-triu( ones(8,'single') ) ), [1 1 1 5] );
                    itd.labelData{2} = repmat( cat(4, single(2)-eye(8,'single'), single(2)-triu( ones(8,'single') ) ), [1 1 1 5] );
                    itd.labelNames = {'L1';'L2'};
                otherwise
                    error(['Invalid sel: ' num2str(sel)]);
            end
            itd.id = 1:numel(itd.set);
        end
        
    end
    
end

        
