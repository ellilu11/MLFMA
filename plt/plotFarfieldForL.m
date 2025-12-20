clear; clc;

maxlvl = 2;

nths = [34,43];
nphs = 2*nths;
nLs = size(nths,2);

flds = cell(nLs,1);
thetas = cell(nLs,1);
phis = cell(nLs,1);

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\ff\";

for i=1:size(nths,2)
    flds{i} = readmatrix(dir+"ff_maxlvl"+string(maxlvl)+"_nth"+string(nths(i))+".txt");
    thetas{i} = readmatrix(dir+"thetas_nth"+string(nths(i))+".txt");
    phis{i} = readmatrix(dir+"phis_nth"+string(nths(i))+".txt");
end

% Select Ls and field components to plot
L1 = 1; L2 = 2;
comp = 2; % E_x = 1, E_y = 2, E_z = 3

%% Reshape fld matrices
fldArr1 = reshape(flds{L1}, nphs(L1), nths(L1), 3);
fldArr2 = reshape(flds{L2}, nphs(L2), nths(L2), 3);

%% Interpolate sampled functions on to regular grid
[ph_1, th_1] = meshgrid(thetas{L1}, phis{L1});
[ph_2, th_2] = meshgrid(thetas{L2}, phis{L2});

nthQ = nths(2)*10;
nphQ = nphs(2)*10;
thetasQ = linspace(0,pi,nthQ);
phisQ = linspace(0,2*pi,nphQ);
[phQ, thQ] = meshgrid(thetasQ, phisQ);

fldInterp1 = interp2(ph_1, th_1, fldArr1(:,:,comp), phQ, thQ, 'linear');
fldInterp2 = interp2(ph_2, th_2, fldArr2(:,:,comp), phQ, thQ, 'linear');

ymax = max([max(fldInterp1), max(fldInterp2)]);

%% Plot as function of theta
close all;
figure(1);
for iph = 1:nphQ
    plot(thetasQ, fldInterp1(iph,:), thetasQ, fldInterp2(iph,:));
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
    plot(phisQ, fldInterp1(:,ith), phisQ, fldInterp2(:,ith));
    xlim([0 2*pi]);
    ylim([0 ymax]);
    xlabel('ph');
    title('th = '+string((ith-1)/nthQ*pi))
    pause(0.01);
end


