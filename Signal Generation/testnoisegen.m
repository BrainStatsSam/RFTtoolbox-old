%A script to test the function noisegen.

Dim = [1,160];
data = noisegen(Dim, 20, 6);

data_mean = mean(data,3);
plot(data_mean);