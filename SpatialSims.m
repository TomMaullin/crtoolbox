% -------------------------------------------------------------------------
% This function takes as inputs:
% -------------------------------------------------------------------------
% - `nSub`: Number of subjects.
% - `simType`: String representing simulation type. e.g. `ramp` for ramp.
% - `nReals`: Number of realizations.
% - `c`: Threshold for mu.
% -------------------------------------------------------------------------
function SpatialSims(nSub, simType, nReals, c)

    % If there is a previous simulation saved, stop as we don't want to
    % overwrite it.
    if exist([simType '.mat'], 'file')
        error('Will not overwrite sim result')
    end
    
    % Tau_n for simulation
    tau = 1/sqrt(nSub);
    
    % Number of bootstraps
    nBoot = 5000;
    
    % Dimension of 2D simulation
    dim = [100 100]; 
    
    % Desired smoothing
    smo = 3;
    
    % When we perform smoothing we extend the edge of the image to avoid
    % edge effects. Following this we remove the edge. The with of the 
    % padding on the edge are given below
    edgeLength = ceil(2*smo/sqrt(2*log(2))); 				 
    
    % Number of voxels
    nVox = prod(dim);
    
    % Adjusted dimensions with edge for smoothing
    adjDim = dim + 2*ceil(edgeLength*smo*ones(1,2));
    
    % Indices for original image with respect to padded image
    xInd = {(ceil(edgeLength*smo)+1):(ceil(edgeLength*smo)+dim(1))};
    yInd = {(ceil(edgeLength*smo)+1):(ceil(edgeLength*smo)+dim(2))};
    
    % Combined indices
    xyInd = cat(2, xInd, yInd);
    
    % This vector stores the result for each realisation on whether AC^+ < AC < AC^ for each level of smoothing (1 if true, 0 if false) 
    obsSuccess_raw_80           = zeros(nRlz, 1); 
    obsSuccess_raw_90           = zeros(nRlz, 1);
    obsSuccess_raw_95           = zeros(nRlz, 1);
    obsSuccess_observed_80           = zeros(nRlz, 1); 
    obsSuccess_observed_90           = zeros(nRlz, 1);
    obsSuccess_observed_95           = zeros(nRlz, 1);
    obsSuccess_raw_80_alternate = zeros(nRlz, 1); 
    obsSuccess_raw_90_alternate = zeros(nRlz, 1);
    obsSuccess_raw_95_alternate = zeros(nRlz, 1);
    obsSuccess_observed_80_alternate = zeros(nRlz, 1); 
    obsSuccess_observed_90_alternate = zeros(nRlz, 1);
    obsSuccess_observed_95_alternate = zeros(nRlz, 1);
    
end


