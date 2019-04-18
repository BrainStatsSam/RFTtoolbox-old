function [ out ] = RFTthresh( cluster_mean )
% RFTthresh calculates the threshold given that the mean number of clusters
% takes a 
%--------------------------------------------------------------------------
% ARGUMENTS
% 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport.

%Need to make a function that computes the expected number of clusters
%above a threshold. And then then write this function so that it acts like
%spm_uc_RF in that it uses Monte Carlo but uses 1 - exp(-Em) where Em is
%the average number of clusters ie the true number! This will allow you to
% verify that the expectedEC value is causing the conservativeness in the
% simulated case.  (And can then extrapolate to the real data case I
% guess?)

%Note that the difference between clusterwise and voxelwise inference
%appears to be the use of the nosko result which makes sense that this
%wouldn't hold as clusters that I have seen are rarely spherical.

%Height is more robust to cluster inference as a small change in the
%threshold could drastically change the size of a cluster but won't affect
%the height much!

end

