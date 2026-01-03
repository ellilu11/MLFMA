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
absFlds = permute(reshape(absFlds, nph, nth, 3), [2 1 3]);

%%
comp = 3;
ymax = max(absFlds(:,:,comp), [], "all");

close all;
figure(1)
for ith = 1:nth
    plot(phis, absFlds(ith,:,comp));
    xlim([0 2*pi]);
    ylim([0 ymax]);
    xlabel('ph');
    title('th = '+string(thetas(ith)));
    pause(0.1);
end

%% 
[theta, phi] = meshgrid(phis, thetas);
X = R*cos(theta).*sin(phi);
Y = R*sin(theta).*sin(phi);
Z = R*cos(phi);

lim = [-rootLeng/2 rootLeng/2];
scale = 2.0;

figure(2);
% quiver3(X, Y, Z, absFlds(:,:,1), absFlds(:,:,2), absFlds(:,:,3))
quiver3InLogScale(X, Y, Z, absFlds(:,:,1), absFlds(:,:,2), absFlds(:,:,3), scale)
xlim(lim); ylim(lim); zlim(lim);
xlabel('x'); ylabel('y'); zlabel('z'); 