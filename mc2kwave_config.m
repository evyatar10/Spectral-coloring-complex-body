function cfg = mc2kwave_config()
% One place for every knob used by the MC->k-Wave bridge.

%% Vessel counts (keep "as before" here)
cfg.vessel_counts.num_superficial_triplets = 30;   % each triplet = 1 artery + 2 veins
cfg.vessel_counts.num_deep_triplets        = 8;

%% Tissue layering & domain (match MCMATLAB defaults)
cfg.layers.water_surface_z_m  = 1.0e-4;   % 0.01 cm
cfg.layers.epidermis_thick_m  = 8.0e-5;   % 0.008 cm (~80 µm)
cfg.layers.tissue_bottom_m    = 1.50e-2;  % 1.5 cm  (MC z-extent)
% bands measured from the water surface (top)
cfg.layers.superficial_band_m = [2.5e-4, 6.0e-4];
cfg.layers.deep_band_m        = [1.0e-3, 2.5e-3];

%% Optical compositions (Jacques model; SAME as your MCMATLAB)
% Values are [S, B, W, F, M].
cfg.optics.epidermis = [0.75, 0.00, 0.75, 0.00, 0.01];
cfg.optics.dermis    = [0.67, 0.002,0.65, 0.00, 0.00];
% Blood base fractions: B=1 fixed, W=0.51, F=0, M=0 (S set per artery/vein)
cfg.optics.blood_WFM = [1.00, 0.51, 0.00, 0.00];

% Artery / vein saturations (defaults if you choose constant model)
cfg.optics.S_artery_default = 0.98;
cfg.optics.S_vein_default   = 0.75;

% Venous sO2 model:
%   'constant'    -> all veins use S_vein_default
%   'layer_dists' -> draw per-vein sO2 from depth-band distributions
cfg.optics.venous_model           = 'layer_dists';
cfg.optics.superficial_vein_mean  = 0.85;  cfg.optics.superficial_vein_std = 0.05;
cfg.optics.deep_vein_mean         = 0.57;  cfg.optics.deep_vein_std        = 0.07;
cfg.optics.artery_mean            = 0.98;  cfg.optics.artery_std           = 0.015;

%% Acoustic & k-Wave settings
cfg.acoust.c0        = 1540;     % m/s
cfg.acoust.rho0      = 1000;     % kg/m^3
cfg.acoust.alpha_db  = 0.75;     % dB/(MHz^y)/cm
cfg.acoust.alpha_pow = 1.5;
cfg.sim.CFL          = 0.30;     % makeTime(...)
cfg.sim.datacast     = 'single';
cfg.sim.make_US      = false;    % set true if/when you add a US sim

%% Sensor layout
cfg.sensor.lateral_coverage = [0.10, 0.90];  % fraction of width at top

%% Randomness (reproducibility)
% Set to a number (e.g., 12345) to make vessel placement & venous sO2 draws repeatable.
% Set to [] to let MATLAB's RNG vary each run.
cfg.random.seed = 12345;
end
