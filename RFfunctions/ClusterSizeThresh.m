function k = ClusterSizeThresh( u, alpha, df, STAT, resel_vec, n )
% ClusterSizeThresh( FWHM ) serves as a function template.
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
if nargin < 2
    alpha = 0.05;
end

dk = 1;
d  = 1000;
k = 0;
while abs(d) > 0
    p         = spm_P_RF(1,k,u,df,STAT,resel_vec,n)
    q         = spm_P_RF(1,k+dk,u,df,STAT,resel_vec,n)
    d         = round((alpha - p)/((q - p)/dk)) %Note d can be negative. %This is just Newton-Raphson!
    k         = k + d
    if isinf(k)
        k=+Inf; 
        return
    end
end

end

