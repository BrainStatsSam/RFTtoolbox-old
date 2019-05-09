function out = Gker( x, sigma2_or_FWHM, use_fwhm )
% GKER( x, sigma2_or_FWHM, use_fwhm ) calculates the Gaussian Kernel given
% data and the variance: sigma2 or FWHM.
%--------------------------------------------------------------------------
% ARGUMENTS
% x
% sigma2_or_FWHM    If FWHM, it is the FWHM in voxels.
% use_fwhm
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport.
if nargin < 3
    use_fwhm = 1;
end

if use_fwhm
    sigma2 = FWHM2sigma(sigma2_or_FWHM);
else
    sigma2 = sigma2_or_FWHM;
end

out = exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2);

end
