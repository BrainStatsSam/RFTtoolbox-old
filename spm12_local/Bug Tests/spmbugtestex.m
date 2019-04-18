%Testing for a bug in spm_conv.m

Sig = SpheroidSignal([256, 256],20,3,0);
Img = spm_conv(Sig,30,30);
surf(Img)

[~, max_index ] = max(Img);

[~, max_index ] = max([1,2;3,4]);
