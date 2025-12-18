clear; clc;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\coeffs\";

coeffs1 = readmatrix(dir+"coeffs.txt");
coeffs2 = readmatrix(dir+"acoeffs.txt");

nvec = 1:length(coeffs1);

%%
comp = 1;

figure(1)
plot(nvec, coeffs1(:,comp), nvec, coeffs2(:,comp));

relErr = abs(coeffs1(:,comp)-coeffs2(:,comp)) ./ coeffs2(:,comp);
figure(2)
semilogy(nvec, relErr)