function [ SigAll ] = gensig_dep( Dim, Smo, Mag, Rad )
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
% %2d Signal
% Sig = gensig_dep(2);
% surf(Sig)
% 
% %3d Signal
% Sig = gensig_dep(3, 6, 2, 10);
% surf(Sig(:,:,50))
%--------------------------------------------------------------------------
% SEE ALSO
% SimDemo

%Set to sum(100*clock) to ensure that this is different each time.
randn('seed',sum(100*clock));   %-Random number generator initializaiton

%DEFAULT VARS
if (nargin < 1)
    Dim = [100,100];
end
if (Dim == 2)
   Dim = [100,100];
end
if (Dim == 3)
    Dim = [91, 109, 91];
end

if (nargin < 2)
    Smo = 6;
end
if (nargin < 3)
    Mag = 2;
end
if (nargin < 4)
    Rad = 10;
end

rimFWHM = 1.7;
%-----------Initialization of Some Variables
nDim    = length(Dim);

nSmo    = length(Smo);       %-Number of smoothing kernel widths

wDim    = Dim + 2*ceil(rimFWHM*max(Smo(:)))*ones(1,nDim);  % Working image dimension

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

% Signal for all smoothings
SigAll = zeros([Dim nSmo]);

% Pre-computations, common to all realisations: makes an ellipsoid of
% signal with no smoothing and no noise added.
Sig = SpheroidSignal(wDim, Rad, Mag, 6); %- Signal
surf(Sig(:,:,50))
pause

for j = 1:nSmo
    %
    % Smooth signal
    %
    if Smo(j)==0
        Sigs      = Sig;
    else
        Sigs      = zeros(wDim);
        if (length(wDim)>2)
            spm_smooth(Sig,Sigs,Smo(:,j));  %Note still in Smo(,) form make need to change this!
            %Also this doesn't calculate Sig.
        else
            [Sigs,~]   = spm_conv(Sig,Smo(j),Smo(j));
            %ss is unnecessary here as we don't need to scale the signal as
            %its not random.
        end
    end
    % Truncate
    if nDim==2
        tSigs    = Sigs(TrnInd{1},TrnInd{2}); %Evaluates the signal on the actual image.
    else
        tSigs    = Sigs(TrnInd{1},TrnInd{2},TrnInd{3});
    end
    % Scale to desired magnitude
    % Q: why is the signal scaled as an input to SpheroidSignal if its scaled
    % again here? Yeah I think unnec to do twice! Probs yes but nec here as
    % convolution has potentially changed the magnitudes.
    tSigs    = Mag*tSigs/max(tSigs(:));
    
    if nDim==2
        SigAll(:,:,j) = tSigs;
    else
        SigAll(:,:,:,j) = tSigs;
    end
    
end

end

