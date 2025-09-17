% demo_run.m  (MC -> k-Wave bridge -> Method A)  with per-stage timers
clear; clc; close all;

addpath(genpath('C:\Users\evyat\MATLAB\MCMATLAB'))
addpath(genpath('C:\Users\evyat\MATLAB\k-Wave'))


% For very-new NVIDIA GPUs with MATLAB: enable CUDA forward-compat
try
    parallel.gpu.enableCUDAForwardCompatibility(true);
    gpuDevice([]); gpuDevice;
catch
    % OK if no GPU / older device; MC code will handle CPU if needed
end

cfg = mc2kwave_config();     % single source of truth

% MC inputs (sizes/depth define the replicated triplets' geometry downstream)
inpMC.wavelengths = 680:10:970;    % nm
inpMC.saturations = [96 75 75];    % [artery, left vein, right vein] (used by MCMATLAB only)
inpMC.Hshape      = 'BOX_1CM_1CM';
inpMC.Diameter    = 0.003;         % m
inpMC.Dist        = 0.004;         % m
inpMC.VLdiam      = 0.0025;        % m
inpMC.VRdiam      = 0.0025;        % m
inpMC.separation  = 0.0005;        % m
inpMC.angle       = 90;            % deg

outdir = fullfile(pwd,'created_data_from_mc');

% ---------- 1) MC (UNCHANGED function) ----------
tMC = tic;
[Hraw, ~] = MCMATLAB(inpMC);
dtMC = toc(tMC);
fprintf('\n[Timing] MCmatlab     : %8.3f s (%.2f min)\n', dtMC, dtMC/60);

% ---------- 2) Bridge: MC -> k-Wave ----------
tKW = tic;
gen_kwave_from_mc(outdir, Hraw, inpMC, cfg);
dtKW = toc(tKW);
fprintf('[Timing] k-Wave bridge: %8.3f s (%.2f min)\n', dtKW, dtKW/60);

% ---------- 3) Your existing analysis (Method A) ----------
tM = tic;
methodA_fixed(fullfile(outdir,'MethodAData.mat'));
dtM = toc(tM);
fprintf('[Timing] Method A      : %8.3f s (%.2f min)\n', dtM, dtM/60);

% ---------- Summary ----------
totalT = dtMC + dtKW + dtM;
fprintf('-----------------------------------------\n');
fprintf('Total pipeline          %8.3f s (%.2f min)\n\n', totalT, totalT/60);
