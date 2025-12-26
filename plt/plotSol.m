clear; clc;

digits = 6;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\sol\";

sol1 = readmatrix(dir+"sol.txt");
% sol1 = readmatrix(dir+"sol_d" + digits +".txt");
sol2 = readmatrix(dir+"solDir.txt");

nvec = 1:length(sol1);

%% 
close all;
for i=1:2
    sol1Sort = sortrows(sol1,i);
    sol2Sort = sortrows(sol2,i);
    figure(2*i-1);
    plot(nvec, sol1Sort(:,i), nvec, sol2Sort(:,i));
    % semilogy(nvec, abs(solSort(:,i)), nvec, abs(solDirSort(:,i)));

    relErr = abs(sol1Sort(:,i)-sol2Sort(:,i)) ./ abs(sol2Sort(:,i));
    figure(2*i);
    semilogy(nvec, relErr)
end

