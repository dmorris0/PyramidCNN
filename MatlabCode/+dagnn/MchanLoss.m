classdef MchanLoss < dagnn.Loss

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
    %multi-channel softmax log loss
    %Assumes inputs{} is size 3
    properties
        type = 'mchan'
        wt = 1     %can weight loss with this
        wtIndices = [1 1]  %indices into global reweight
    end
    
    methods
        %constructor:
        function obj = MchanLoss(varargin)
            obj.load(varargin) ;
        end
        
        function outputs = forward(obj, inputs, params)
            [sz1,sz2,~] = size(inputs{1});
            scale = obj.wt * single(1/(sz1*sz2));  %normalize by number of pixels per image
            outputs{1} = vl_nnmchansoftmaxloss( inputs{:}, [], obj.wtIndices ) * scale;
            n = obj.numAveraged ;
            m = n + size(inputs{1},4) ;
            obj.average = (n * obj.average + gather(outputs{1})) / m ;
            obj.numAveraged = m ;
        end
        
        function [derInputs, derParams] = backward(obj, inputs, params, derOutputs)
            [sz1,sz2,~] = size(inputs{1});
            scale = obj.wt * single(1/(sz1*sz2));  %normalize by number of pixels per image
            derInputs{1} = vl_nnmchansoftmaxloss(inputs{:}, derOutputs{1}, obj.wtIndices) * scale ;
            derInputs{2} = [] ;
            derInputs{3} = [] ;
            derParams = {} ;
        end
        
        function outputSizes = getOutputSizes(obj, inputSizes, paramSizes)
            outputSizes{1} = [1 1 1 inputSizes{1}(4)] ;
        end
        
        function rfs = getReceptiveFields(obj)
            % the receptive field depends on the dimension of the variables
            % which is not known until the network is run
            rfs(1,1).size = [NaN NaN] ;
            rfs(1,1).stride = [NaN NaN] ;
            rfs(1,1).offset = [NaN NaN] ;
            rfs(2,1) = rfs(1,1) ;
            rfs(3,1) = rfs(1,1) ;
        end
        
    end
    
end
