function gen_kwave_network(save_dir)
% gen_kwave_network(save_dir)
% Photoacoustic sim with:
%   - Single, laterally-uniform surface illumination (no dual beams)
%   - Realistic many-vessel vascular bed (superficial + deep plexuses)
%   - Epidermal melanin layer (avascular), dermis/hypodermis with background Hb
%   - Same output schema as your gen_kwave_fixed -> Method A runs unchanged
%
% Output: writes <save_dir>/MethodAData.mat with RF, PAI, US, parameters, energy...

if nargin < 1
    save_dir = fullfile(pwd,'created_data_network');
end
if ~exist(save_dir,'dir'), mkdir(save_dir); end

%% Grid & acoustic medium (kept close to your defaults where possible)
Nx = 128; Nz = 128;
dx = 0.1e-3; dz = 0.1e-3;                 % 0.1 mm pixels (100 µm) to keep runtime reasonable
[x_grid, z_grid] = meshgrid((0:Nx-1)*dx, (0:Nz-1)*dz);

medium.sound_speed = 1540;                % m/s
medium.density     = 1000;                % kg/m^3

kgrid = kWaveGrid(Nz, dz, Nx, dx);
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed, 0.3);
Fs = 1/dt;

%% Wavelength set
wavelength = 680:10:970;                  % nm
W = numel(wavelength);

%% Layers & masks (gel + slab; epidermis avascular, dermis vascular)
gel             = 0.5e-3;                 % 0.5 mm gel/stand-off
slab_thickness  = 8e-3;                   % 8 mm tissue
mask_tissue     = (z_grid >= gel) & (z_grid <= gel + slab_thickness);

% Simple epidermis thickness ~0.15 mm inside the tissue (avascular)
% (papillary plexus lies ~0.3–0.35 mm from surface on average)
epi_thickness = 0.15e-3;                  % 150 µm
mask_epi      = mask_tissue & (z_grid <= gel + epi_thickness);
mask_dermis   = mask_tissue & ~mask_epi;

%% Build a many-vessel vascular bed (superficial + deep plexuses)
opts = struct();
opts.superficial_z_mm = [0.25 0.60];      % mm from gel → ~subpapillary band
opts.deep_z_mm        = [1.0 2.5];        % mm from gel → deep plexus band
opts.num_superficial  = 40;               % small vessels
opts.num_deep         = 15;                % larger vessels
opts.superficial_r_mm = [0.10 0.20];      % radii in mm
opts.deep_r_mm        = [0.30 0.80];      % radii in mm
opts.artery_so2       = 0.95; %0.99
opts.vein_so2         = 0.7;
opts.artery_fraction  = 0.5;              % ~half arteries, half veins
[mask_vess, so2_map]  = local_make_vascular_bed(x_grid, z_grid, gel, opts);

%% Optical absorption & single-source fluence
% Uses standard Hb spectra + melanin/water/fat references (Prahl/OMLC, Jacques)
[mu_a, phi] = local_buildMuAandFluence_network( ...
    wavelength, x_grid, z_grid, mask_epi, mask_dermis, mask_vess, so2_map, gel);

%% Initial pressure p0
Gamma      = 1;           % Grüneisen (unitless)
F0_mJcm2   = 10;          % incident radiant exposure
mu_a_SI    = mu_a * 1000; % mm^-1 -> m^-1
F0_SI      = F0_mJcm2 * 10;       % mJ/cm^2 -> J/m^2
F_SI       = F0_SI * phi;         % same spatial dims as mu_a
p0_all     = Gamma * (mu_a_SI .* F_SI);

%% Linear sensor array (top 20–80% width)
min_sensor_position = 0.1; max_sensor_position = 0.9;
sensor.mask = false(Nz, Nx);
col1 = max(1, round(min_sensor_position*Nx)); col2 = min(Nx, round(max_sensor_position*Nx));
sensor.mask(1, col1:col2) = true;

%% k-Wave forward (per wavelength)
numSensors = sum(sensor.mask(:));
numTime    = numel(kgrid.t_array);
sensor_data = zeros(numSensors, numTime, W, 'single');

for k = 1:W
    source.p0 = p0_all(:,:,k);
    sd = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PlotSim', false, 'DataCast','single');
    sensor_data(:,:,k) = sd;
end

% [time x sensors x wavelengths]
RF = cell(1,1);
RF{1,1} = permute(sensor_data, [2 1 3]);

%% Simple PA back-projection image per wavelength (for your overview plots)
cols   = find(sensor.mask(1,:));  x_sens = (cols-1) * dx;
x_vec  = (0:Nx-1) * dx;  z_vec = (0:Nz-1) * dz;
PAI_recon = zeros(Nz, Nx, W, 'single');
for k = 1:W
    rf = squeeze(sensor_data(:,:,k));     % Ns x Nt
    PAI_recon(:,:,k) = local_pa_backproj(rf, Fs, medium.sound_speed, x_vec, z_vec, x_sens);
end
PAI = { PAI_recon };

%% Ultrasound B-mode (reuses your helper, but with many vessels)
fc_us     = 5e6;      % Hz
numCycles = 2;
[USimg, ~] = local_us_kwave_plane(kgrid, medium, mask_tissue, mask_vess, sensor.mask, dx, dz, fc_us, numCycles);
US = { repmat(single(USimg), [1 1 W]) };

%% Parameters & save
% % parameters.PApixelmmHeight = medium.sound_speed * dt * 1e3;  % mm/sample
% parameters.PApixelmmHeight = dz * 1e3;        % mm per image pixel (correct for plots)
% parameters.PApixelmmWidth  = dx * 1e3;                       % mm/pixel
% parameters.PAheight        = size(RF{1,1},1);
% parameters.PAwidth         = size(RF{1,1},2);
% parameters.PAdepthOffset   = 0;

parameters.PArfmmHeight  = medium.sound_speed * dt * 1e3;   % RF depth step (mm/sample)  ← for RF plots
parameters.PAimgmmHeight = dz * 1e3;                        % IMAGE depth step (mm/pixel)← for image/sO2 plots
parameters.PApixelmmHeight = parameters.PArfmmHeight;       % keep old name for backward-compat (RF)

parameters.PApixelmmWidth  = dx * 1e3;                      % lateral pixel (mm)
parameters.PAdepthOffset   = 0;                             % (optional) use gel/epi offsets if you want

parameters.gel_mm          = gel * 1e3;
parameters.epidermis_mm    = epi_thickness * 1e3;
parameters.tissue_mm       = slab_thickness * 1e3;


energy = ones(1, W, 'double');

save(fullfile(save_dir,'MethodAData.mat'), 'wavelength','Fs','PAI','RF','US','parameters','energy','-v7.3');
fprintf('Saved %s\n', fullfile(save_dir,'MethodAData.mat'));

% Store the spectra struct so your analysis can reuse exactly the same basis
A = getHemoglobinMu(wavelength); 
save(fullfile(save_dir,'mu_invivo.mat'), 'A');  % same filename as your code expects, but contains A
end

% ---------- LOCAL HELPERS (kept compatible with your style) ----------

function PAI = local_pa_backproj(rf, Fs, c, x_grid, z_grid, x_sens)
[Ns, Nt] = size(rf);
t = (0:Nt-1)/Fs; Nx = numel(x_grid); Nz = numel(z_grid);
apo = hann(Ns).'; PAI = zeros(Nz, Nx);
for ix = 1:Nx
    x = x_grid(ix);
    for iz = 1:Nz
        z = z_grid(iz);
        r   = hypot(x - x_sens, z);
        tau = r/c;
        s = zeros(1,Ns);
        for sIdx = 1:Ns
            s(sIdx) = interp1(t, rf(sIdx,:), tau(sIdx), 'linear', 0);
        end
        PAI(iz,ix) = sum(apo .* s);
    end
end
PAI = max(PAI,0);
end

function [USimg, us_rf] = local_us_kwave_plane(kgrid, medium_base, mask_tissue, mask_vess, sensor_mask, dx, dz, fc, numCycles)
% Same as your helper, just takes many-vessel mask.
[Nz, Nx] = size(mask_tissue);
c0   = single(medium_base.sound_speed);
rho0 = single(medium_base.density);

c   = c0  * ones(Nz, Nx, 'single');
rho = rho0* ones(Nz, Nx, 'single');

% weak random heterogeneity to generate backscatter
if exist('imgaussfilt','file')
    speck = imgaussfilt(randn(Nz, Nx, 'single'), 2);
else
    H = fspecial('gaussian',[9 9],2);
    speck = imfilter(randn(Nz, Nx, 'single'), H, 'replicate');
end
c(mask_tissue)   = c(mask_tissue)   + 15*speck(mask_tissue);
rho(mask_tissue) = rho(mask_tissue) .* (1 + 0.01*speck(mask_tissue));

% vessel acoustic contrast + boundary
c(mask_vess)   = c0*0.985;
rho(mask_vess) = rho0*0.985;
ring = bwperim(mask_vess,8);
c(ring)   = c0*1.02; rho(ring) = rho0*1.02;

medium.sound_speed  = c;
medium.density      = rho;
medium.alpha_coeff  = 0.75;     % dB/(MHz^y)/cm
medium.alpha_power  = 1.5;

% TX plane-wave from full array, RX on same
cols = find(sensor_mask(1,:)); Ns = numel(cols);
src_mask = false(size(sensor_mask)); src_mask(1, cols) = true;
Fs_sim    = 1 / kgrid.dt;
input_sig = toneBurst(Fs_sim, fc, numCycles, 'Envelope','Gaussian');
source.p_mask = src_mask;
source.p      = repmat(input_sig, [sum(src_mask(:)), 1]);
sensor.mask   = sensor_mask;

args = {'PMLInside', false, 'PlotSim', false, 'DataCast','single'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, args{:});
us_rf = sensor_data;

% DAS two-way
x_sens = (cols-1) * dx; x_grid = (0:Nx-1)*dx; z_grid = (0:Nz-1)*dz;
USrf_img = local_us_das(us_rf, Fs_sim, double(c0), x_grid, z_grid, x_sens);

% envelope + log compression
env = local_envelope_img(USrf_img); env = env / max(env(:) + eps);
USimg_db = 20*log10(env + eps);                   % [-60,0] dB
USimg = 200 * (USimg_db - (-60)) / 60;            % ~[0..200]
USimg = single(max(min(USimg, 200), 0));
end

function USimg = local_us_das(rf, Fs, c, x_grid, z_grid, x_sens)
[Ns, Nt] = size(rf); t = (0:Nt-1)/Fs; Nx = numel(x_grid); Nz = numel(z_grid);
apo = hann(Ns).'; USimg = zeros(Nz, Nx);
for ix = 1:Nx
    x = x_grid(ix);
    for iz = 1:Nz
        z = z_grid(iz);
        r   = hypot(x - x_sens, z);
        tau = 2*r/c;
        s = zeros(1,Ns);
        for sIdx = 1:Ns
            s(sIdx) = interp1(t, rf(sIdx,:), tau(sIdx), 'linear', 0);
        end
        USimg(iz,ix) = sum(apo .* s);
    end
end
USimg = max(USimg,0);
end

function env = local_envelope_img(rf_img)
[Nz, Nx] = size(rf_img);
Nt = Nz; H = zeros(Nt,1);
if mod(Nt,2)==0, H(1)=1; H(Nt/2+1)=1; H(2:Nt/2)=2; else, H(1)=1; H(2:(Nt+1)/2)=2; end
F = fft(rf_img, [], 1);
env = abs(ifft(F .* H, [], 1));
end
