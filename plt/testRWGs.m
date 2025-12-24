clear; clc;

% Import RWGs
dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\config\n1\";

vertices = readmatrix(dir+"vertices.txt");
faces = readmatrix(dir+"faces.txt");
rwgs = readmatrix(dir+"rwgs.txt");

% Specify plane wave
th = pi/4; ph = pi/4;
[x,y,z] = sph2cart(1.0, th, ph); kvec = [x,y,z];

