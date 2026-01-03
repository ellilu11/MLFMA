clear; clc;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\";

nsrcs = 1000;
srcs = readmatrix(dir+"config\dipole\sphere_n"+nsrcs+".txt");

sols = readmatrix(dir+"out\sol\solDir.txt");
assert(length(sols) == nsrcs);

rootLeng = 10;
R = 0.9*rootLeng/2;

%% 
srcx = srcs(:,1);
srcy = srcs(:,2);
srcz = srcs(:,3);
assert(max(sqrt(srcx.^2 + srcy.^2 + srcz.^2)-R) < 1.0E-4); 

vals = sqrt(sols(:,1).^2 + sols(:,2).^2); % complex modulus of sols

npts = 100;
[theta, phi] = meshgrid(linspace(0,2*pi,2*npts), linspace(0,pi,npts));
X = R.*cos(theta).*sin(phi);
Y = R.*sin(theta).*sin(phi);
Z = R.*cos(phi);

F_val = scatteredInterpolant(srcx, srcy, srcz, vals, 'linear');
Vq = F_val(X, Y, Z);

close all;
figure(1);
s = surf(X, Y, Z, Vq);

s.FaceColor = 'interp'; 
s.EdgeColor = 'none';

%%
% figure(2)
% scatter3(srcx, srcy, srcz)