clear; % clc;

digits = 6;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\sol\";

sol1 = readmatrix(dir+"rvec.txt");
sol2 = readmatrix(dir+"rvec_interp.txt");
solDir = readmatrix(dir+"rvecDir.txt");
% sol1 = readmatrix(dir+"rvec_nq7.txt");
% sol2 = readmatrix(dir+"rvec_nq7_interp.txt");
% solDir = readmatrix(dir+"rvecDir_nq7.txt");
% sol1 = readmatrix(dir+"curr_nq7.txt");
% sol2 = readmatrix(dir+"curr_nq7_interp.txt");
% solDir = readmatrix(dir+"currDir_nq7.txt");

nvec = 1:length(sol1);

%% 
close all;
for i=1:2
    sol1Sort = sortrows(sol1,i);
    sol2Sort = sortrows(sol2,i);
    solDirSort = sortrows(solDir,i);
    figure(2*i-1);
    % plot(nvec, sol1Sort(:,i), nvec, solDirSort(:,i));
    semilogy(nvec, abs(sol1Sort(:,i)), nvec, abs(solDirSort(:,i)));

    figure(2*i);

    err1 = abs(sol1Sort(:,i)-solDirSort(:,i));
    err2 = abs(sol2Sort(:,i)-solDirSort(:,i));
    semilogy(nvec, err1, nvec, err2);

    % relErr1 = abs(sol1Sort(:,i)-solDirSort(:,i)) ./ abs(solDirSort(:,i));
    % relErr2 = abs(sol2Sort(:,i)-solDirSort(:,i)) ./ abs(solDirSort(:,i));
    % semilogy(nvec, relErr1, nvec, relErr2);
end

%%
mean(err1)
mean(err2)

% mean(relErr1)
% mean(relErr2)
