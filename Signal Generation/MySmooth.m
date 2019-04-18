function sImg = MySmooth(Img,smo)
% MYSMOOTH smooths the image Img using the methods from spm_smooth. It also 
% normalizes the resulting image so that the highest value of the signal is 
% 1.
%--------------------------------------------------------------------------
% ARGUMENTS
% Img   An array of 1s and 0s. 1s indicate where we would like the signal
%       to be. Can be a 2d or a 3d array.
% smo   A 1x2 or 1x3 vector that is the smoothness in terms of the FWHM, 
%       ie smo = [FWHM_x, FWHM_y, FWHM_z]: the FWHM in each of the x, y and
%       z directions.
%--------------------------------------------------------------------------
% OUTPUT
% sImg  An array of the same dimensions of Img that has been smoothed
%       according to the degree of smoothness specified by smo.
%--------------------------------------------------------------------------
% EXAMPLES
% 2d version
% sImg = MySmooth(double(MkRadImg([256, 256],[128.5,128.5]) <= 20), 10);
% surf(sImg);
%
% 3d version
% sImg = MySmooth(double(MkRadImg([256, 256, 256],[128.5,128.5, 128.5]) <= 20), 10);
% surf(sImg(:,:,128));
% surf(sImg(:,:,109));
%--------------------------------------------------------------------------
% SEE ALSO
% spm_conv, spm_smooth

Dim     = size(Img); %Calculate the dimensions of the Image.
nDim    = length(Dim); %Calculate the number of dimensions.

%If only 1 smoothing parameter given smooth the same in the x and y directions.
if length(smo) == 1
    smo = repmat(smo, 1, nDim); 
end

if any(smo)
    if nDim == 2
        sImg = spm_conv(double(Img), smo(1), smo(2));
    elseif nDim == 3
        sImg = 0*Img; %Create an empty image of the right dimensions.
        spm_smooth(double(Img),sImg,smo); 
    else
        error('smo must be a 1x2 or a 1x3 vector of positive numbers.');
    end
    sImg = sImg/max(sImg(:)); 
else
    %Note need double here to convert from logical to a vector.
    sImg = Img/max(double(Img(:)));
end

%Note that you divide by the maximum above so that you can control the
%magnitude of the signal externally to this function. So mag in the parent
%function SpheroidSignal is the value that the highest signal takes.

%This seems a bit arbitrary.
%sImg(sImg(:)<0.05) = 0;

return