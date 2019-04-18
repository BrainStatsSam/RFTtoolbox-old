function Mean = SimDemo(nSubj,nRlz,SvNm)
% Simulation Demo: It calls SpheroidSignal.m to generate a spheroid of 
% signal.
%
% nSubj = Number of subjects.
% SvNm  = File name to be saved.
% nRlz  = Number of toy runs, is the number of times you simulate noise for
% each subject.
%----------------------------------------------
%

%Set to sum(100*clock) to ensure that this is different each time.
randn('seed',sum(100*clock));   %-Random number generator initializaiton

%------------Starting Up initialization
if (nargin<1)
  nSubj  = 60;  % Number of subjects
end
if (nargin<2)
  nRlz = 5000;
end  
if (nargin<3)
  SvNm  = 'Normsim';  % Save name
end
% if exist([SvNm '.mat'], 'file')
%   error('Will not overwrite sim result')
% end

%------------Define parameters

Dim    = [256, 256];  	  %-Spatial extent of signal 

Smo = [1.5, 3];

%% If you want to do 3D:
% Dim   = [112 112 16];  	  %-Spatial extent of signal   
% 
% Smo   = [1.5 3 6];               %-Smoothness in terms of FHWM

Rad   = 20; 			  %-Equatorial radius of spheroid

Mag   = 3;			  %-Magnitude of signal

rimFWHM = 3;       % The width of padding around the image (xFWHM) to be 
                   % truncated to reduce the edge effect.

%-----------Initialization of Some Variables
nDim    = length(Dim);

nSmo    = length(Smo);       %-Number of smoothing kernel widths

%Smo     = ones(nDim,1)*Smo;     %-Smoothness (isotropic) (vector of ones)

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

% Index for truncation
V       = prod(Dim);                          %-Image volume (voxels)

% All noise (unsmoothed) for all subjects, allowing an outer smoothing loop
RawNoise = zeros([wDim nSubj]); %This creates an array of zeros.

% Signal+Noise for all subjects, for a given smoothing level (reused in
% side smoothing loop)
Data   = zeros([Dim nSubj]);

% Signal for all smoothings
SigAll = zeros([Dim nSmo]);

% Pre-computations, common to all realisations: makes an ellipsoid of
% signal with no smoothing and no noise added.
Sig = SpheroidSignal(wDim, Rad, Mag, 0); %- Signal


for j = 1:nSmo
  % 
  % Smooth signal 
  %
  if Smo(j)==0
    Sigs      = Sig;
  else
    Sigs      = zeros(wDim);
    if (length(wDim)>2)
      ss       = spm_smooth(Sig,Sigs,Smo(:,j)');  %Note still in Smo(,) form make need to change this!
      %Also this doesn't calculate Sig.
    else
      [Sigs,ss]   = spm_conv(Sig,Smo(j),Smo(j));
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
  % again here? Yeah I think unnec to do twice!
  tSigs    = Mag*tSigs/max(tSigs(:));
  
  if nDim==2
    SigAll(:,:,j) = tSigs; 
  else
    SigAll(:,:,:,j) = tSigs; 
  end

end      

%
% Realisations loop
%
for t=1:nRlz
  for j= 1:nSmo

    fprintf('Realization %5d FWHM %d mm\n',t,Smo(1,j));
    
    for i=1:nSubj
      %
      % Generate random realizations (once for all smoothings, ie only j=1.)
      %
      if j==1
        if nDim==2   %-Independent N(0,1) noise that will be added to the signal.
            RawNoise(:,:,i)   = randn(wDim); 
        else
            RawNoise(:,:,:,i)   = randn(wDim);
        end
      end
      
      %
      % Smooth noise  
      %
      if Smo(1,j)==0
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
            tt       = spm_smooth(RawNoise(:,:,:,i),Noises,Smo(:,j)');
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
      
      if nDim==2
        Data(:,:,i) = tImgs;
      else
        Data(:,:,:,i) = tImgs;
      end

    end %========== Loop i (subjects)
    
    Mean = mean(Data,nDim+1); %This returns the mean over a number of subjects.

    % Do some other stuff...


  end  %-------- Loop j (Different levels of smoothing)
end %------ Loop t (nRlz)

%eval(['save ' SvNm ' nSubj nRlz Dim wDim Smo Rad Mag rimFWHM ']); % whatever else

