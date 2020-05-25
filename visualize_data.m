clc
clear all
close all



u = load('step60000_centraldiff_u_data.txt');
v = load('step60000_centraldiff_v_data.txt');
phi = load('step60000_centraldiff_phi_data.txt');
figure
quiver(u,v)
figure
surf(phi)