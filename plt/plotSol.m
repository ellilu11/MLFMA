clear; clc;

digits = 6;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\sol\";

sol1 = readmatrix(dir+"curr_nq7.txt");
sol2 = readmatrix(dir+"currDir_nq7.txt");
% sol1 = readmatrix(dir+"rvec_nq7.txt");
% sol2 = readmatrix(dir+"rvecDir_nq7.txt");
% sol1 = readmatrix(dir+"dip\curr.txt");
% sol2 = readmatrix(dir+"dip\currDir.txt");

nvec = 1:length(sol1);

%% 
close all;
for i=1:2
    sol1Sort = sortrows(sol1,i);
    sol2Sort = sortrows(sol2,i);
    figure(2*i-1);
    plot(nvec, sol1Sort(:,i), nvec, sol2Sort(:,i));
    % semilogy(nvec, abs(sol1Sort(:,i)), nvec, abs(sol2Sort(:,i)));

    relErr = abs(sol1Sort(:,i)-sol2Sort(:,i)) ./ abs(sol2Sort(:,i));
    figure(2*i);
    semilogy(nvec, relErr)
end

