function [inter_smoothed_data, ss_vec ] = inter_conv1D(data, increment, FWHM, normalize)
% INTER_CONV1D performs convolution on a set of data points where the
% interior values are calculated as well.
%--------------------------------------------------------------------------
% ARGUMENTS
% increment     the constant amount that each voxel is partioned into.
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% [inter_smoothed_data, ss_vec ] = inter_conv1D([1:4], 1, 2)
% [inter_smoothed_data, ss_vec ] = inter_conv1D([1:4], 0.5, 2)
% [inter_smoothed_data, ss_vec ] = inter_conv1D([1:4], 0.1, 2)
%--------------------------------------------------------------------------
% NOTES
% Potentially should change that so that (given length n) instead of 
% generating data on [1,n] it generates it on [1/2, n+1/2] which might make
% more sense in the context of the brain.
% AUTHOR: Sam Davenport.
if nargin < 4
    normalize = 1;
end

%Ensure that the data is a 1D row vector.
s = size(data);
if length(s) > 2
    error('Need a 1D array')
elseif s(2) < s(1)
    s = fliplr(s);
    data = data';
end
if s(1) ~= 1
    error('Need a 1D array')
end

kernel_param = FWHM/sqrt(8*log(2)) + eps; %Why: so you don't divide by zero.
range_of_conv = min(round(6*kernel_param), s(2));
% range_of_conv = round(6*kernel_param);
x = -range_of_conv:range_of_conv;

inter_smoothed_data = 1:increment:length(data);
inter_smoothed_data = [inter_smoothed_data; zeros(1, length(inter_smoothed_data))];

set_of_increms = 0:increment:(1-increment);
nincrems = length(set_of_increms);

ss_vec = zeros(1, length(inter_smoothed_data));
%The case I = 1 must be covered separately include the final data point as
%well.
kernel    = exp(-(x-set_of_increms(1)).^2/(2*kernel_param^2))/sqrt(2*pi*kernel_param^2);
smoothed_data = conv(data, kernel); %This has length length(data) + length(kernel) - 1
inter_smoothed_data( 2, 1:nincrems:(nincrems*length(data))) = smoothed_data((range_of_conv + 1):end - range_of_conv); %this cuts off the first bit of incomplete overlap
ss_vec(1:nincrems:(nincrems*length(data))) = sum(kernel.^2);

for I = 2:length(set_of_increms)
    kernel    = fliplr(exp(-(x-set_of_increms(I)).^2/(2*kernel_param^2))/(sqrt(2*pi)));
    %     kernel    = kernel/sum(kernel); %You shouldn't scale this as that
    %     will set things to be out of whack. (technical term I'll have you
    %     know). This is the K we're using see Fabian's RFT in fMRI.
    smoothed_data = conv(data, kernel); %This has length length(data) + length(kernel) - 1
    inter_smoothed_data( 2, I:nincrems:(nincrems*(length(data)-1))) = smoothed_data((range_of_conv+1):end - range_of_conv - 1); %need to change the range here
    % have range_of_conv on one side and on the other so need to get rid of the
    % first and last lot (plus one extra as the interior elements have one
    % element less than length(data) does.
    ss_vec(I:nincrems:(nincrems*(length(data)-1))) = sum(kernel.^2);
end

if normalize == 1
    inter_smoothed_data(2, :) = inter_smoothed_data(2, :)./sqrt(ss_vec);
end

end

