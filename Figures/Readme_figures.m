rftboxloc = '/data/fireback/davenpor/davenpor/Toolboxes/RFTtoolbox/';

clf
pos_vector = [0,550,1000,600];
set(gcf, 'position', pos_vector)

noise = noisegen(160,20,6);
plot(mean(noise,2), 'linewidth', 2)
xlabel('Voxels')
ylabel('Mean over 20 subjects')
title('1D noise generation')

export_fig([rftboxloc, 'Figures/readme_1Dreal.png'], '-transparent')


clf
pos_vector = [0,550,2000,1000];
set(0,'defaultAxesFontSize', 20);
set(gcf, 'position', pos_vector)

Dim = [160,160];
noise = noisegen(Dim, 20, 6);
noise_mean = mean(noise,3);
surf(noise_mean)
title('2D noise generation')
xlabel('x')
ylabel('y')

export_fig([rftboxloc, 'Figures/readme_2Dreal.png'], '-transparent')

clf
pos_vector = [0,550,2000,1000];
set(0,'defaultAxesFontSize', 20);
set(gcf, 'position', pos_vector)

title('2D Signal Generation')
xlabel('x')
ylabel('y')

Sig = gensig([1.3,2], 3, [10,20], [100,150], {[40,30], [70,120]});
surf(Sig)

export_fig([rftboxloc, 'Figures/readme_signal.png'], '-transparent')
