clear; clc;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\";

sol = readmatrix(dir+"valsInterped.txt");
solDir = readmatrix(dir+"valsDirect.txt");

nvec = 1:length(sol);

%% 
close all;
for i=1:2
    figure(i);
    plot(nvec, sol(:,i), nvec, solDir(:,i));

    % relErr = abs(sol(:,i)-solDir(:,i)) ./ abs(solDir(:,i));
    % figure(2*i);
    % semilogy(nvec, relErr)
end
