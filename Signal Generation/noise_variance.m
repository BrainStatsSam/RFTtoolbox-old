%% 1D
nvox = 1000;
nsubj = 20;
data = zeros(nsubj, nvox);
for I = 1:nsubj
    noise = randn(nvox,1);
    [smoothed_noise, ss] = spm_conv(noise, 6);
    data(I, :) = smoothed_noise./sqrt(ss);
    plot(data(I,:), 'linewidth', 2)
end

mean(std(data,0 ,1))

%% 2D
Dim = [100,100];
nsubj = 20;
nvox = prod(Dim);
data = zeros(nsubj, nvox);
for I = 1:nsubj
    noise = randn(Dim(1), Dim(2));
    [smoothed_noise, ss] = spm_conv(noise, 6);
    data(I, :) = smoothed_noise(:)./sqrt(ss);
end

mean(std(data,0 ,1))

%%
% noise = noisegen(160, nsubj, 6);
% s = std(noise, 0,2);
% 
% %%
% noise = noisegen([50,50], 50, 6, 1);
% s = std(noise, 0,2);

