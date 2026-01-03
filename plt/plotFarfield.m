clear; clc;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\ff\";

flds = readmatrix(dir+"ff_nq7.txt");
thetas = readmatrix(dir+"thetas.txt");
phis = readmatrix(dir+"phis.txt");

nth = length(thetas); nph = length(phis);
assert(length(flds) == nth*nph);

rootLeng = 20.0;
R = 0.9*rootLeng/2;

%% 
absFlds = abs(flds(:,1:2:5) + 1i*flds(:,2:2:6));
absFlds = reshape(absFlds, nth, nph, 3);

[theta, phi] = meshgrid(phis, thetas);
X = R*cos(theta).*sin(phi);
Y = R*sin(theta).*sin(phi);
Z = R*cos(phi);

%% 
lim = [-rootLeng/2 rootLeng/2];

figure(1);
quiver3(X, Y, Z, absFlds(:,:,1), absFlds(:,:,2), absFlds(:,:,3))
xlim(lim); ylim(lim); zlim(lim);
xlabel('x'); ylabel('y'); zlabel('z'); 