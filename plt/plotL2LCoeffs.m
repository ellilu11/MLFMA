clear; clc;

nths = [11,16];
nphs = 2*nths;
nLs = size(nths,2);

coeffs = cell(nLs,1);
thetas = cell(nLs,1);
phis = cell(nLs,1);

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\coeffs\";

for i=1:size(nths,2)
    coeffs{i} = readmatrix(dir+"lcoeffs_nth"+string(nths(i))+".txt");
    thetas{i} = readmatrix(dir+"thetas_nth"+string(nths(i))+".txt");
    phis{i} = readmatrix(dir+"phis_nth"+string(nths(i))+".txt");
end

% Select Ls and field components to plot
L1 = 1; L2 = 2;
comp = 2;

%% Reshape fld matrices
coeffArr1 = reshape(coeffs{L1}, nphs(L1), nths(L1), 6);
coeffArr2 = reshape(coeffs{L2}, nphs(L2), nths(L2), 6);

%% Interpolate sampled functions on to regular grid
[ph_1, th_1] = meshgrid(thetas{L1}, phis{L1});
[ph_2, th_2] = meshgrid(thetas{L2}, phis{L2});

nthQ = nths(2)*10;
nphQ = nphs(2)*10;
thetasQ = linspace(0,pi,nthQ);
phisQ = linspace(0,2*pi,nphQ);
[phQ, thQ] = meshgrid(thetasQ, phisQ);

coeffInterp1 = interp2(ph_1, th_1, coeffArr1(:,:,comp), phQ, thQ, 'linear');
coeffInterp2 = interp2(ph_2, th_2, coeffArr2(:,:,comp), phQ, thQ, 'linear');

ymax = max([max(coeffInterp1), max(coeffInterp2)]);

%% Plot as function of theta
close all;
figure(1);
for iph = 1:nphQ
    plot(thetasQ, coeffInterp1(iph,:), thetasQ, coeffInterp2(iph,:));
    xlim([0 pi]);
    ylim([0 ymax]);
    xlabel('th');
    title('ph = '+string((iph-1)/nphQ*2*pi))
    pause(0.01);
end

%% Plot as function of phi
close all;
figure(2);
for ith = 1:nthQ
    plot(phisQ, coeffInterp1(:,ith), phisQ, coeffInterp2(:,ith));
    xlim([0 2*pi]);
    ylim([0 ymax]);
    xlabel('ph');
    title('th = '+string((ith-1)/nthQ*pi))
    pause(0.01);
end


