clear; clc;

digits = 6;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\sol\";

sol = readmatrix(dir+"sol_d" + digits +".txt");
% sol = readmatrix(dir+"solDir_recip.txt");
solDir = readmatrix(dir+"solDir.txt");

nvec = 1:length(sol);

%% 
close all;
for i=1:2
    solSort = sortrows(sol,i);
    solDirSort = sortrows(solDir,i);
    figure(2*i-1);
    % plot(nvec, solSort(:,i), nvec, solDirSort(:,i));
    semilogy(nvec, abs(solSort(:,i)), nvec, abs(solDirSort(:,i)));

    relErr = abs(solSort(:,i)-solDirSort(:,i)) ./ abs(solDirSort(:,i));
    figure(2*i);
    semilogy(nvec, relErr)
end

