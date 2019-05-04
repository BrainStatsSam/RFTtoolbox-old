function Sig = gensig( Mag, Rad, Smo, Dim, centre_locs )
% GENSIG( Dim, Smo, Mag, Rad, rimFWHM ) generates smoothed signal.
%--------------------------------------------------------------------------
% ARGUMENTS
% nSubj     is the number of subjects which corresponds to the number of
%           separate images to generate.
% Dim       A 1 by 2 or 1 by 3 vector of the image dimensions,
%           Dim = [256,256,256] corresponds to a 256*256*256 image.
% Rad       A positive number that is the equatorial radius of Spheroid.
% Mag       the magnitude of the signal.
% Smo       a real number that equals FWHM_x = FWHM_y = FWHM_z ie we
%           assume that the smoothing is equal in all directions. This can
%           also be a vector in which case it will specify different levels
%           of smoothing. ATM smoothes the signal with a given smoothness.
%           Something to change probably or maybe not.
% rimFWHM   the width of padding around the image to be truncated to reduce
%           the edge effect.
%--------------------------------------------------------------------------
% OUTPUT
% data      an array of size Dim(1) by Dim(2) by nSubj where the third
%           index runs over the number of subject: for each subject giving
%           an image, points of which are identified by the first two
%           indices.
%--------------------------------------------------------------------------
% EXAMPLES
% %3D signal
% Sig = gensig(2, 10, 6);
% surf(Sig(:,:,50))
%--------------------------------------------------------------------------
% SEE ALSO
% SimDemo

%Set to sum(100*clock) to ensure that this is different each time.
randn('seed',sum(100*clock));   %-Random number generator initializaiton

%DEFAULT VARS

if nargin < 1
    Mag = 2;
end
if nargin < 2
    Rad = 10;
end
if nargin < 3
    Smo = 6;
end
if nargin < 4
    Dim = [91, 109, 91];
end
if nargin < 5
    centre_locs = {Dim/2 + 1/2};
end
    
npeaks = length(centre_locs);
if length(Mag) == 1
    Mag = repmat(Mag, 1, npeaks);
elseif length(Mag) ~= npeaks
    error('The number of peaks in Mag is not the same as in centre_peaks')
end

if length(Rad) == 1
    Rad = repmat(Rad, 1, npeaks);
elseif length(Rad) ~= npeaks
    error('The number of peaks in Rad is not the same as in centre_peaks')
end

if length(Smo) == 1
    Smo = repmat(Smo, 1, npeaks);
elseif length(Smo) ~= npeaks
    error('The number of peaks in Smo is not the same as in centre_peaks')
end

rimFWHM = 1.7;
%-----------Initialization of Some Variables
nDim    = length(Dim);

boundary2add = ceil(rimFWHM*max(Smo(:)))*ones(1,nDim);
wDim    = Dim + 2*boundary2add;  % Working image dimension

for peak = 1:npeaks
    centre_locs{peak} = centre_locs{peak} + boundary2add; % Need to do this to ensure the centre locs are in the right place even after boundary stuff.
end

%Below describes the actual bits of the image other than the extra
%truncation stuff that we have added on!
Trunc_x = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(1))};
Trunc_y = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(2))};

if nDim==2
    %Concatenates Trunc_x and Trunc_y into one array. Why is this
    %necessary? 
    TrnInd = cat(2, Trunc_x, Trunc_y); %Note cat(2,A,B) == [A,B]
else
    Trunc_z = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(3))};
    TrnInd  = cat(2, Trunc_x, Trunc_y, Trunc_z);
end

Sig = zeros(wDim);
for peak = 1:npeaks
    Sig = Sig + SpheroidSignal(wDim, Rad(peak), Mag(peak), Smo(peak), centre_locs{peak}); %- Signal Should smooth here really!!
end

if nDim == 2
    Sig = Sig(TrnInd{1},TrnInd{2});
elseif nDim == 3
    Sig = Sig(TrnInd{1},TrnInd{2}, TrnInd{3});
end

end

