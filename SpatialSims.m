% -------------------------------------------------------------------------
% This function takes as inputs:
% -------------------------------------------------------------------------
% - `nSub`: Number of subjects.
% - `simType`: String representing simulation type. e.g. `ramp` for ramp.
% - `nReals`: Number of realizations.
% - `c`: a for mu.
% -------------------------------------------------------------------------
function SpatialSims(nSub, simType, nReals, c)

    tic
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
    
    % MARKER: I THINK THIS CAN BE REMOVED BUT NEED TO INVESTIGATE
    stdblk  = prod(dim([1 2])/2);

    % When we perform smoothing we extend the edge of the image to avoid
    % edge effects. Following this we remove the edge. The with of the 
    % padding on the edge are given below
    edgeLength = ceil(2*smo/sqrt(2*log(2))); 	 
    
    % Adjusted dimensions with edge for smoothing
    adjDim = dim + 2*ceil(edgeLength*smo*ones(1,2));			 
    
    % Number of voxels
    nVox = prod(dim);
    
    % Indices for original image with respect to padded image
    xInd = {(ceil(edgeLength*smo)+1):(ceil(edgeLength*smo)+dim(1))};
    yInd = {(ceil(edgeLength*smo)+1):(ceil(edgeLength*smo)+dim(2))};
    
    % Combined indices
    xyInd = cat(2, xInd, yInd);
    
    % Observed booleans for \hat{A}_c^+ < A_c < \hat{A}_c^- (i.e. using the
    % true boundary) assessed via gridpoints
    obsSuccess_trueBdry_80 = zeros(nReals, 1); % MARKER: THIS CAN BE MADE 1 LINE INSTEAD OF 3 
    obsSuccess_trueBdry_90 = zeros(nReals, 1);
    obsSuccess_trueBdry_95 = zeros(nReals, 1);
    
    % Observed booleans for \hat{A}_c^+ < \hat{A}_c < \hat{A}_c^- (i.e. 
    % using the estimated boundary) assessed via gridpoints
    obsSuccess_estBdry_80 = zeros(nReals, 1); % MARKER: THIS CAN BE MADE 1 LINE INSTEAD OF 3
    obsSuccess_estBdry_90 = zeros(nReals, 1);
    obsSuccess_estBdry_95 = zeros(nReals, 1);
    
    % Observed booleans for \hat{A}_c^+ < A_c < \hat{A}_c^- (i.e. using the
    % true boundary) assessed via interpolation
    obsSuccess_trueBdry_80_interpol = zeros(nReals, 1); % MARKER: THIS CAN BE MADE 1 LINE INSTEAD OF 3
    obsSuccess_trueBdry_90_interpol = zeros(nReals, 1);
    obsSuccess_trueBdry_95_interpol = zeros(nReals, 1);
    
    % Observed booleans for \hat{A}_c^+ < \hat{A}_c < \hat{A}_c^- (i.e. 
    % using the estimated boundary) assessed via interpolation
    obsSuccess_estBdry_80_interpol = zeros(nReals, 1); % MARKER: THIS CAN BE MADE 1 LINE INSTEAD OF 3
    obsSuccess_estBdry_90_interpol = zeros(nReals, 1);
    obsSuccess_estBdry_95_interpol = zeros(nReals, 1);

    % This vector stores the a value for each run, obtained from a boostrap
    % on the true boundary
    a_trueBdry_80 = zeros(nReals, 1); % MARKER: THIS CAN BE MADE 1 LINE INSTEAD OF 3
    a_trueBdry_90 = zeros(nReals, 1);
    a_trueBdry_95 = zeros(nReals, 1);

    % This vector stores the a value for each run, obtained from a boostrap
    % on the estimated boundary
    a_estBdry_80 = zeros(nReals, 1); % MARKER: THIS CAN BE MADE 1 LINE INSTEAD OF 3
    a_estBdry_90 = zeros(nReals, 1);
    a_estBdry_95 = zeros(nReals, 1);

    % This is where we will save the boostrap results for the supremum of
    % G
    supG_trueBdry = zeros(nBoot,1);
    supG_estBdry = zeros(nBoot,1);
    
    % Create mu
    mu = create_signal(dim, simType, 1);
    
    % Work out the true A_c
    Ac = mu >= c;

    % Boundary edges of true A_c
    Ac_bdry_edges = getBdryparams(mu, c);
    toc 
    disp('Marker 1')

    % =====================================================================
    % Loop through realisations
    % =====================================================================
    for t=1:nReals
        
        tic
        % Print a dot
        fprintf('.');
        
        % Initialize the estBdry data
        observed_data = zeros([dim nSub]);

        % =================================================================
        % Loop through subjects (data generation)
        % =================================================================
        % For each subject generate an image
        for i=1:nSub
            
            % Generate noise
            epsilon = create_noise(adjDim, 'homo', 1, smo, xyInd); % MARKER: Prelim testing shows smoothing can be done on 3d array slice by slice
            
            % Generate 
            img = mu + epsilon; 
            
            % Save estBdry data
            observed_data(:,:,i) = img;

        end 

        % Take the mean of the data
        muHat = mean(observed_data,3);

        % Take the standard deviation of the date
        est_std = reshape(biasmystd(reshape(observed_data,...
                                                [nVox nSub]),...
                                         stdblk),... % MARKER: This reshaping looks dubious - Investigate
                               dim);

        toc
        disp('Marker 2')
                           
        % Work out the estimated \hat{A}_c
        AcHat = muHat >= c;
        
        % Boundary edges of estimated \hat{A}_c
        AcHat_bdry_edges = getBdryparams(muHat, c);

        % Residuals (i.e. data minus mean)
        resid = bsxfun(@minus,observed_data,muHat);
        
        % Standardise the residuals (MARKER: Would be good to go through
        % this operation)
        resid = spdiags(1./reshape(est_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*reshape(resid,[prod(dim) nSub]); 

        % =================================================================
        % Loop through Bootstrap instances
        % =================================================================
        % Implementing the Multiplier Boostrap to obtain confidence intervals
        tic
        for k=1:nBoot 

            % Rademacher variables/signflips for this bootstrap
            signflips = randi(2,[nSub,1])*2-3;

            % Obtain bootstrap residuals
            resid_bootstrap = resid*spdiags(signflips, 0, nSub, nSub);
            resid_bootstrap = reshape(resid_bootstrap, [dim nSub]);
            
            % Obtain G bootstrap field
            G_bootstrap_field = sum(resid_bootstrap, 3)/sqrt(nSub);
            
            % Re-standardizing by bootstrap standard deviation
            boot_std = std(resid_bootstrap, 0, 3);
            G_bootstrap_field = G_bootstrap_field./boot_std; % MARKER: This isn't standard practice but Alex said they found it worked well for small sample. Left in for now.

            % Calculating the maximum over the weighted interpolated true boundary edges
            trueBdry_values = getBdryvalues(G_bootstrap_field, Ac_bdry_edges);
            supG_trueBdry(k) = max(abs(trueBdry_values)); 

            % Calculating the maximum over the weighted interpolated estimated boundary edges
            estBdry_values = getBdryvalues(G_bootstrap_field, AcHat_bdry_edges);
            supG_estBdry(k) = max(abs(estBdry_values));   

        end
        toc 
        disp('Marker 3')
        
        % =================================================================
        % Set computations using Ac boundary: 
        % -----------------------------------------------------------------
        % Here Achatminus\Ac and Ac\Achatplus are computed
        % =================================================================
        
        tic
        % Get the a values along the true boundary for each percentile 
        % MARKER: Unnecessary repetition
        a_trueBdry_80 = prctile(supG_trueBdry, 80);
        a_trueBdry_90 = prctile(supG_trueBdry, 90);
        a_trueBdry_95 = prctile(supG_trueBdry, 95);
        
        % Get \hat{A}_c^+ and \hat{A}_c^- for p=0.8 
        AcHatminus_trueBdry_80 = muHat >= c - a_trueBdry_80*tau*sigma;
        AcHatplus_trueBdry_80 = muHat >= c + a_trueBdry_80*tau*sigma;
        
        % Get \hat{A}_c^{+}\A_c and A_c\(\hat{A}_c^{-})
        AcHatplus_setminus_Ac_trueBdry_80 = AcHatplus_trueBdry_80.*(1 - Ac);
        Ac_setminus_AcHatminus_trueBdry_80 = Ac.*(1 - AcHatminus_trueBdry_80);
        
        % Get \hat{A}_c^+ and \hat{A}_c^- for p=0.9 
        AcHatminus_trueBdry_90 = muHat >= c - a_trueBdry_90*tau*sigma; % MARKER: Unnecessary repetition
        AcHatplus_trueBdry_90 = muHat >= c + a_trueBdry_90*tau*sigma;
        
        % Get \hat{A}_c^{+}\A_c and A_c\(\hat{A}_c^{-})
        AcHatplus_setminus_Ac_trueBdry_90 = AcHatplus_trueBdry_90.*(1 - Ac);
        Ac_setminus_AcHatminus_trueBdry_90 = Ac.*(1 - AcHatminus_trueBdry_90);
        
        % Get \hat{A}_c^+ and \hat{A}_c^- for p=0.95 
        AcHatminus_trueBdry_95 = muHat >= c - a_trueBdry_95*tau*sigma; % MARKER: Unnecessary repetition
        AcHatplus_trueBdry_95 = muHat >= c + a_trueBdry_95*tau*sigma;
        
        % Get \hat{A}_c^{+}\A_c and A_c\(\hat{A}_c^{-})
        AcHatplus_setminus_Ac_trueBdry_95 = AcHatplus_trueBdry_95.*(1 - Ac);
        Ac_setminus_AcHatminus_trueBdry_95 = Ac.*(1 - AcHatminus_trueBdry_95);
        toc
        disp('Marker 4')

        % =================================================================
        % Set computations using \hat{Ac} boundary: 
        % -----------------------------------------------------------------
        % Here Achatminus\Ac and Ac\Achatplus are computed
        % =================================================================
        
        % Get the a values along the estimated boundary for each percentile 
        % MARKER: Unnecessary repetition
        a_estBdry_80 = prctile(supG_estBdry, 80);
        a_estBdry_90 = prctile(supG_estBdry, 90);
        a_estBdry_95 = prctile(supG_estBdry, 95);

        % Get \hat{A}_c^+ and \hat{A}_c^- for p=0.8 
        AcHatminus_estBdry_80 = muHat >= c - a_estdry_80*tau*sigma;
        AcHatplus_estBdry_80 = muHat >= c + a_estBdry_80*tau*sigma;
        
        % Get \hat{A}_c^{+}\A_c and A_c\(\hat{A}_c^{-}) (Note: We still use
        % true A_c as this is for coverage testing)
        AcHatplus_setminus_Ac_estBdry_80 = AcHatplus_estBdry_80.*(1 - Ac);
        Ac_setminus_AcHatminus_estBdry_80 = Ac.*(1 - AcHatminus_estBdry_80);
        
        % Get \hat{A}_c^+ and \hat{A}_c^- for p=0.9 
        AcHatminus_estBdry_90 = muHat >= c - a_estBdry_90*tau*sigma; % MARKER: Unnecessary repetition
        AcHatplus_estBdry_90 = muHat >= c + a_estBdry_90*tau*sigma;
         
        % Get \hat{A}_c^{+}\A_c and A_c\(\hat{A}_c^{-})
        AcHatplus_setminus_Ac_estBdry_90 = AcHatplus_estBdry_90.*(1 - Ac);
        Ac_setminus_AcHatminus_estBdry_90 = Ac.*(1 - AcHatminus_estBdry_90);

        % Get \hat{A}_c^+ and \hat{A}_c^- for p=0.95 
        AcHatminus_estBdry_95 = muHat >= c - a_estBdry_95*tau*sigma; % MARKER: Unnecessary repetition
        AcHatplus_estBdry_95 = muHat >= c + a_estBdry_95*tau*sigma;
         
        % Get \hat{A}_c^{+}\A_c and A_c\(\hat{A}_c^{-})
        AcHatplus_setminus_Ac_estBdry_95 = AcHatplus_estBdry_95.*(1 - Ac);
        Ac_setminus_AcHatminus_estBdry_95 = Ac.*(1 - AcHatminus_estBdry_95);
        
        % =================================================================
        % Calculate values along boundarys based on supG derived using Ac
        % =================================================================
        
%         % Get the threshold field for c-a\tau\sigma(s)
%         muHat_lowThresh_trueBdry_80 = c - a_trueBdry_80*tau*sigma; % MARKER: Unnecessary repitition
%         muHat_lowThresh_trueBdry_90 = c - a_trueBdry_90*tau*sigma;
%         muHat_lowThresh_trueBdry_95 = c - a_trueBdry_95*tau*sigma;
%         
%         % Get the threshold field for c+a\tau\sigma(s)
%         muHat_uppThresh_trueBdry_80 = c + a_trueBdry_80*tau*sigma;
%         muHat_uppThresh_trueBdry_90 = c + a_trueBdry_90*tau*sigma;
%         muHat_uppThresh_trueBdry_95 = c + a_trueBdry_95*tau*sigma;

        
%         lower_condition_80_boundary_values = getBdryvalues(muHat_lowThresh_trueBdry_80, cohen_d_boundary_edges);
%         upper_condition_80_boundary_values = getBdryvalues(muHat_uppThresh_trueBdry_80, cohen_d_boundary_edges);
% 
%         lower_condition_90_boundary_values = getBdryvalues(muHat_lowThresh_trueBdry_90, cohen_d_boundary_edges);
%         upper_condition_90_boundary_values = getBdryvalues(muHat_uppThresh_trueBdry_90, cohen_d_boundary_edges);
% 
%         lower_condition_95_boundary_values = getBdryvalues(muHat_lowThresh_trueBdry_95, cohen_d_boundary_edges);
%         upper_condition_95_boundary_values = getBdryvalues(muHat_uppThresh_trueBdry_95, cohen_d_boundary_edges);
% 
%         estBdry_cohen_d_true_boundary_values = getBdryvalues(estBdry_cohen_d, cohen_d_boundary_edges);
% 
%         lower_condition_80_success = estBdry_cohen_d_true_boundary_values < lower_condition_80_boundary_values;
%         upper_condition_80_success = estBdry_cohen_d_true_boundary_values >= upper_condition_80_boundary_values;
% 
%         lower_condition_90_success = estBdry_cohen_d_true_boundary_values < lower_condition_90_boundary_values;
%         upper_condition_90_success = estBdry_cohen_d_true_boundary_values >= upper_condition_90_boundary_values;
% 
%         lower_condition_95_success = estBdry_cohen_d_true_boundary_values < lower_condition_95_boundary_values;
%         upper_condition_95_success = estBdry_cohen_d_true_boundary_values >= upper_condition_95_boundary_values;
% 
% 
%         lower_condition_80_estBdry_boundary_values = getBdryvalues(lower_condition_80_estBdry, cohen_d_boundary_edges);
%         upper_condition_80_estBdry_boundary_values = getBdryvalues(upper_condition_80_estBdry, cohen_d_boundary_edges);
% 
%         lower_condition_90_estBdry_boundary_values = getBdryvalues(lower_condition_90_estBdry, cohen_d_boundary_edges);
%         upper_condition_90_estBdry_boundary_values = getBdryvalues(upper_condition_90_estBdry, cohen_d_boundary_edges);
% 
%         lower_condition_95_estBdry_boundary_values = getBdryvalues(lower_condition_95_estBdry, cohen_d_boundary_edges);
%         upper_condition_95_estBdry_boundary_values = getBdryvalues(upper_condition_95_estBdry, cohen_d_boundary_edges);
% 
%         lower_condition_80_estBdry_success = estBdry_cohen_d_true_boundary_values < lower_condition_80_estBdry_boundary_values;
%         upper_condition_80_estBdry_success = estBdry_cohen_d_true_boundary_values >= upper_condition_80_estBdry_boundary_values;
% 
%         lower_condition_90_estBdry_success = estBdry_cohen_d_true_boundary_values < lower_condition_90_estBdry_boundary_values;
%         upper_condition_90_estBdry_success = estBdry_cohen_d_true_boundary_values >= upper_condition_90_estBdry_boundary_values;
% 
%         lower_condition_95_estBdry_success = estBdry_cohen_d_true_boundary_values < lower_condition_95_estBdry_boundary_values;
%         upper_condition_95_estBdry_success = estBdry_cohen_d_true_boundary_values >= upper_condition_95_estBdry_boundary_values;

        % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
        % binarized sets for residuals on the true boundary in mult. bootstrap
        if sum(AcHatplus_setminus_Ac_trueBdry_80(:))+sum(Ac_setminus_AcHatminus_trueBdry_80(:))==0
            obsSuccess_trueBdry_80(t) = 1;
            fprintf('trueBdry nominal 80 success! \n');
        else 
            obsSuccess_trueBdry_80(t) = 0; 
            fprintf('trueBdry nominal 80 failure! \n');
        end

        if sum(AcHatplus_setminus_Ac_trueBdry_90(:))+sum(Ac_setminus_AcHatminus_trueBdry_90(:))==0
            obsSuccess_trueBdry_90(t) = 1;
            fprintf('trueBdry nominal 90 success! \n');
        else 
            obsSuccess_trueBdry_90(t) = 0; 
            fprintf('trueBdry nominal 90 failure! \n');
        end

        if sum(AcHatplus_setminus_Ac_trueBdry_95(:))+sum(Ac_setminus_AcHatminus_trueBdry_95(:))==0
            obsSuccess_trueBdry_95(t) = 1;
            fprintf('trueBdry nominal 95 success! \n');
        else 
            obsSuccess_trueBdry_95(t) = 0; 
            fprintf('trueBdry nominal 95 failure! \n');
        end

        % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
        % binarized sets for residuals on the estBdry boundary in mult. bootstrap
        if sum(AcHatplus_setminus_Ac_estBdry_80(:))+sum(Ac_setminus_AcHatminus_estBdry_80(:))==0
            obsSuccess_estBdry_80(t) = 1;
            fprintf('estBdry nominal 80 success! \n');
        else 
            obsSuccess_estBdry_80(t) = 0; 
            fprintf('estBdry nominal 80 failure! \n');
        end

        if sum(AcHatplus_setminus_Ac_estBdry_90(:))+sum(Ac_setminus_AcHatminus_estBdry_90(:))==0
            obsSuccess_estBdry_90(t) = 1;
            fprintf('estBdry nominal 90 success! \n');
        else 
            obsSuccess_estBdry_90(t) = 0; 
            fprintf('estBdry nominal 90 failure! \n');
        end

        if sum(AcHatplus_setminus_Ac_estBdry_95(:))+sum(Ac_setminus_AcHatminus_estBdry_95(:))==0
            obsSuccess_estBdry_95(t) = 1;
            fprintf('estBdry nominal 95 success! \n');
        else 
            obsSuccess_estBdry_95(t) = 0; 
            fprintf('estBdry nominal 95 failure! \n');
        end

%         % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
%         % binarized sets as well as the linear interpolated boundary method for
%         % residuals taken along the true boundary
%         if sum(AcHatplus_setminus_Ac_trueBdry_80(:))+sum(Ac_setminus_AcHatminus_trueBdry_80(:)+sum(upper_condition_80_success)+sum(lower_condition_80_success))==0
%             subset_success_vector_trueBdry_80_alternate(t) = 1;
%             fprintf('trueBdry nominal 80 alternate true boundary success! \n');
%         else 
%             subset_success_vector_trueBdry_80_alternate(t) = 0; 
%             fprintf('trueBdry nominal 80 alternate true boundary failure! \n');
%         end 
% 
%         if sum(AcHatplus_setminus_Ac_trueBdry_90(:))+sum(Ac_setminus_AcHatminus_trueBdry_90(:)+sum(upper_condition_90_success)+sum(lower_condition_90_success))==0
%             subset_success_vector_trueBdry_90_alternate(t) = 1; 
%             fprintf('trueBdry nominal 90 alternate true boundary success! \n');
%         else 
%             subset_success_vector_trueBdry_90_alternate(t) = 0; 
%             fprintf('trueBdry nominal 90 alternate true boundary failure! \n');
%         end 
% 
%         if sum(AcHatplus_setminus_Ac_trueBdry_95(:))+sum(Ac_setminus_AcHatminus_trueBdry_95(:)+sum(upper_condition_95_success)+sum(lower_condition_95_success))==0
%             subset_success_vector_trueBdry_95_alternate(t) = 1; 
%             fprintf('trueBdry nominal 95 alternate true boundary success! \n');
%         else 
%             subset_success_vector_trueBdry_95_alternate(t) = 0; 
%             fprintf('trueBdry nominal 95 alternate true boundary failure! \n');
%         end 
% 
%         % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
%         % binarized sets as well as the linear interpolated boundary method for
%         % residuals taken along the estBdry boundary
%         if sum(AcHatplus_setminus_Ac_estBdry_80(:))+sum(Ac_setminus_AcHatminus_estBdry_80(:)+sum(upper_condition_80_estBdry_success)+sum(lower_condition_80_estBdry_success))==0
%             subset_success_vector_estBdry_80_alternate(t) = 1;
%             fprintf('estBdry nominal 80 alternate true boundary success! \n');
%         else 
%             subset_success_vector_estBdry_80_alternate(t) = 0; 
%             fprintf('estBdry nominal 80 alternate true boundary failure! \n');
%         end 
% 
%         if sum(AcHatplus_setminus_Ac_estBdry_90(:))+sum(Ac_setminus_AcHatminus_estBdry_90(:)+sum(upper_condition_90_estBdry_success)+sum(lower_condition_90_estBdry_success))==0
%             subset_success_vector_estBdry_90_alternate(t) = 1; 
%             fprintf('estBdry nominal 90 alternate true boundary success! \n');
%         else 
%             subset_success_vector_estBdry_90_alternate(t) = 0; 
%             fprintf('estBdry nominal 90 alternate true boundary failure! \n');
%         end 
% 
%         if sum(AcHatplus_setminus_Ac_estBdry_95(:))+sum(Ac_setminus_AcHatminus_estBdry_95(:)+sum(upper_condition_95_estBdry_success)+sum(lower_condition_95_estBdry_success))==0
%             subset_success_vector_estBdry_95_alternate(t) = 1; 
%             fprintf('estBdry nominal 95 alternate true boundary success! \n');
%         else 
%             subset_success_vector_estBdry_95_alternate(t) = 0; 
%             fprintf('estBdry nominal 95 alternate true boundary failure! \n');
%         end 

    end

    

end

% addpath('C:\Users\user\Downloads\spm12\')
% addpath('C:\Users\user\Documents\Confidence_Sets_Manuscript\custom_functions')