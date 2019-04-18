function [ data ] = datagen( Dim, nSubj, Smo, Mag, Rad, rimFWHM )
% DATAGEN( Dim, nSubj, Smo, Mag, Rad, rimFWHM ) generates 2d smoothed 
% images with signal and correlated gaussian noise that comes from smoothed
% N(0,1) random variables.
%--------------------------------------------------------------------------
% ARGUMENTS
% nSubj     is the number of subjects which corresponds to the number of
%           separate images to generate.
% Dim       Gives different options that correspond to the dimensions.
%           5: 2D noise. 6: 3D noise.
% Rad       A positive number that is the equatorial radius of Spheroid.
% Mag       the magnitude of the signal.
% Smo       a real number that equals FWHM_x = FWHM_y = FWHM_z ie we
%           assume that the smoothing is equal in all directions. This can
%           also be a vector in which case it will specify different levels
%           of smoothing.
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
% data = datagen(2);
% surf(mean(data, 3));
%
% data = datagen(3, 20);
% mean_est = mean(data, 4);
% surf(mean_est(:,:,50));
% surf(mean_est(:,:,30));
%
% 2D random noise:
% data = datagen(5);
% surf(mean(data,3))
%
% % 3D noise with signal
% data = datagen(stdsize, 20, 6);
% m_data = mean(data, 4);
% surf(m_data(:,:,50))
%--------------------------------------------------------------------------
% SEE ALSO
% spm_conv, spm_smooth, SpheroidSignal, 

%Set to sum(100*clock) to ensure that this is different each time.
randn('seed',sum(100*clock));   %-Random number generator initializaiton

%DEFAULT VARS
reshape = 0;
if (nargin < 1)
    Dim = [100,100];
end
if (Dim == 2)
   Dim = [100,100];
end
if (Dim == 3)
    Dim = [91, 109, 91];
end
if (Dim == 4)
    reshape = 1;
    Dim = [91, 109, 91];
end
use_sig = 1;
if (Dim == 5)
    %Generate random noise images. 2D.
    use_sig = 0;
    Dim = [100, 100];
end
if (Dim == 6)
    %Generate random noise images. 3D.
    use_sig = 0;
    Dim = [91, 109, 91];
end

if (nargin < 2)
    nSubj  = 20;  % Number of subjects
end
if (nargin < 3)
    Smo = 6;
end
if (nargin < 4)
    Mag = 2;
end
if (nargin < 5)
    Rad = 10;
end
if (nargin < 6)
    rimFWHM = 3;
end

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

% All noise (unsmoothed) for all subjects, allowing an outer smoothing loop
RawNoise = zeros([wDim nSubj]); %This creates an array of zeros.

% Signal+Noise for all subjects, for a given smoothing level (reused in
% side smoothing loop)
data   = zeros([Dim nSubj]);
if reshape
    data = zeros(nSubj, prod(Dim));
end

% Signal for all smoothings
SigAll = zeros([Dim nSmo]);

% Pre-computations, common to all realisations: makes an ellipsoid of
% signal with no smoothing and no noise added.

if use_sig == 1
    Sig = SpheroidSignal(wDim, Rad, Mag, 0); %- Signal
    for j = 1:nSmo
        %
        % Smooth signal
        %
        if Smo(j)==0
            Sigs      = Sig;
        else
            %This bit doesn't like it when Sig is 0!!! annoyingly. Need to go over
            %this.
            if (length(Dim) == 3)
                Sigs      = zeros(wDim);
                spm_smooth(Sig,Sigs,Smo(j)); %A C++ function so it changes Sigs
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
else
    SigAll = zeros([Dim, nSmo]);
end
%
% Realisations loop
%

for j= 1:nSmo
    for i=1:nSubj
        %
        % Generate random realizations (once for all smoothings, ie only j=1.)
        %
        if j==1
            if nDim == 2   %-Independent N(0,1) noise that will be added to the signal.
                RawNoise(:,:,i)   = randn(wDim);
            else
                RawNoise(:,:,:,i)   = randn(wDim);
            end
        end
        
        %
        % Smooth noise
         %
        if Smo(j)==0
            if (nDim==2)
                Noises    = RawNoise(:,:,i);
            else
                Noises    = RawNoise(:,:,:,i);
            end
        else
            Noises    = zeros(wDim);
            if (nDim==2)
                [Noises,tt] = spm_conv(RawNoise(:,:,i),Smo(j),Smo(j));
            else
                tt       = spm_smooth(RawNoise(:,:,:,i),Noises,Smo(:,j));
                %C++ code, for saying that RawNoise is taken in and Noises
                %is the smooth noise that is the output.
            end
            Noises    = Noises/sqrt(tt); %Done to standardize.
        end
        
        %
        % Truncate to avoid edge effects, create signal+noise
        %
        if nDim==2
            tNoises    = Noises(TrnInd{1},TrnInd{2});
            tImgs      = SigAll(:,:,j) + tNoises;
        else
            tNoises    = Noises(TrnInd{1},TrnInd{2},TrnInd{3});
            tImgs      = SigAll(:,:,:,j) + tNoises;
        end
        
        % Store smoothed truncated signal image in Sigdatamat and signal +
        % noise image in Data
        
        if reshape
            data(i,:) = tImgs(:);
        else
            if nDim==2
            elseif nDim==2
                data(:,:,i) = tImgs;
            elseif  nDim==3
                data(:,:,:,i) = tImgs;
            end
        end
        
    end %========== Loop i (subjects)
    
    %Mean = mean(Data,nDim+1); %This returns the mean over a number of subjects.
    
end  %-------- Loop j (Different levels of smoothing)


%eval(['save ' SvNm ' nSubj nRlz Dim wDim Smo Rad Mag rimFWHM ']); % whatever else
end

