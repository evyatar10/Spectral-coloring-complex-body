% demo_run_network.m  (optional helper)
clear; close all; clc;
addpath(genpath('C:\Users\evyatarr\Documents\MATLAB\k-Wave'));

outdir = fullfile(pwd,'created_data_network');
gen_kwave_network(outdir);                                  % NEW: many-vessel, single-source sim
methodA_fixed(fullfile(outdir,'MethodAData.mat'));          % unchanged Method A
