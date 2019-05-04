function rfs = genRF( nreal, df, FWHM, Dim, asvector )
% genRF( nreal, df, FWHM, Dim, asvector ) returns an nsubj by prod(Dim) 
% set of random fields with degrees of freedom df and a given FWHM. 
%--------------------------------------------------------------------------
% ARGUMENTS
% nreal:the number of realizations
% df    If df == 1, a GRF is generated. If df = [1, n] a t field with
%       n degrees of freedom is generated. If df = [m, n] an F field 
%       with m and n degrees of freedom is generated.
% FWHM  The smoothness.
% Dim   The dimensions of each field. Default is Dim = [91,109,91].
%       asvector determines whether an nreal by prod(Dim) or a nreal by Dim 
%       array is returned. Default takes as_vector = 1.
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% rfs = genRF(10, [1,5], 5)
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport.
if nargin < 4
    Dim = [91,109,91];
end
if nargin < 5
    asvector = 1;
end
if asvector ~= 1
    error('Need to check through the t generation as haven''t set an option for this yet')
end

if asvector == 1
    rfs = zeros(nreal, prod(Dim));
else
    rfs = zeros(nreal, Dim);
end

if isequal(df, 1)
    if asvector == 1;
        rfs = noisegen( Dim, nreal, FWHM, 3 );
    else
        rfs = noisegen( Dim, nreal, FWHM, 1 );
    end
else
    if length(df) ~= 2 
        error('Misoptioning')
    end
end

if df(1) == 1 && length(df) > 1
    for real = 1:nreal
        grfs = noisegen( Dim, (df(2)+1), FWHM, 3 );
        [~,~,~,rfs(real, :)] = meanmos(grfs); %Generates the t-statistic.
    end
elseif df(1) > 1
    error('Haven''t designed generation of F fields yet, ie where the df in the numerator is not 1!')
    for real = 1:nreal
        grfs = noisegen( Dim, (df(2)+1), FWHM, 3 );
        [~,~,~,rfs(real, :)] = meanmos(grfs); %Generates the t-statistic.
    end
end

end

