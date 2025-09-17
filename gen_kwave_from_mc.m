function gen_kwave_from_mc(save_dir, Hraw, inpMC, cfg)
% Build k-Wave inputs from MCMATLAB output using a single config (cfg).
% No constants duplicated — everything comes from cfg or inpMC.

if nargin < 4, cfg = mc2kwave_config(); end
if nargin < 1 || isempty(save_dir), save_dir = fullfile(pwd,'created_data_from_mc'); end
if ~exist(save_dir,'dir'), mkdir(save_dir); end
assert(~isempty(Hraw),'Hraw from MCMATLAB is required');
assert(isstruct(inpMC),'inpMC is required');

% RNG (for reproducible vessel placement & venous sO2, if desired)
if ~isempty(cfg.random.seed), rng(cfg.random.seed); end

wavelength = inpMC.wavelengths(:)';  L = numel(wavelength);

% --- reshape MC output to [Nz x Nx x L]
if ndims(Hraw)==4
    Hxz = squeeze(mean(Hraw,2));   % [Nx x Nz x L]
elseif ndims(Hraw)==3
    Hxz = Hraw;                    % [Nx x Nz x L]
else
    error('Hraw must be (Nx x Ny x Nz x L) or (Nx x Nz x L)');
end
[Nx, Nz, Lh] = size(Hxz);
assert(Lh==L,'Hraw wavelengths mismatch inpMC.wavelengths');
H_mc = permute(Hxz,[2 1 3]);       % [Nz x Nx x L]

% --- MC physical domain (fixed in your MCMATLAB)
Lx_m = 0.02;  Lz_m = 0.015;
dx = Lx_m / Nx;  dz = Lz_m / Nz;
[x_grid, z_grid] = meshgrid((0:Nx-1)*dx, (0:Nz-1)*dz);

% --- k-Wave grid & medium
kgrid = kWaveGrid(Nz, dz, Nx, dx);
medium.sound_speed = cfg.acoust.c0;
medium.density     = cfg.acoust.rho0;
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed, cfg.sim.CFL);
Fs = 1/dt;
medium.alpha_coeff = cfg.acoust.alpha_db;
medium.alpha_power = cfg.acoust.alpha_pow;

% --- masks from cfg layers
z0   = cfg.layers.water_surface_z_m;
zepi = cfg.layers.epidermis_thick_m;
zbot = cfg.layers.tissue_bottom_m;
mask_water  = (z_grid < z0);
mask_epi    = (z_grid >= z0) & (z_grid < z0 + zepi);
mask_dermis = (z_grid >= z0 + zepi) & (z_grid <= zbot);

% --- sanity: ensure MC tissue z-extent >= configured tissue_bottom
if abs(Lz_m - zbot) > 5*dz
    warning('cfg.layers.tissue_bottom_m (%.3g m) != MC z-extent (%.3g m). Using MC extent.', ...
        zbot, Lz_m);
end

% --- optical compositions (MATCH your MCMATLAB defaults)
[S_e,B_e,W_e,F_e,M_e] = unpack(cfg.optics.epidermis);
[S_d,B_d,W_d,F_d,M_d] = unpack(cfg.optics.dermis);
B_blood = 1.0; W_blood = cfg.optics.blood_WFM(2); F_blood = 0; M_blood = 0;

mua_epi = cm_to_m(calc_mua(wavelength, S_e, B_e, W_e, F_e, M_e));  % 1xL
mua_der = cm_to_m(calc_mua(wavelength, S_d, B_d, W_d, F_d, M_d));  % 1xL

% --- background μa (no vessels)
mu_a_bg = zeros(Nz,Nx,L,'single');
for k = 1:L
    mu = zeros(Nz,Nx,'single');
    mu(mask_water)  = cm_to_m(0.00036);  % water μa (cm^-1 -> m^-1)
    mu(mask_epi)    = single(mua_epi(k));
    mu(mask_dermis) = single(mua_der(k));
    mu_a_bg(:,:,k)  = mu;
end

% --- derive vessel radii/offsets from YOUR MCMATLAB inputs (sizes/depth from MC)
rA  = inpMC.Diameter/2;
rLV = inpMC.VLdiam/2;
rRV = inpMC.VRdiam/2;
zcA = inpMC.Dist;
Roff= rA + max(rLV,rRV) + inpMC.separation;
dxA = Roff * sind(inpMC.angle);
dzA = Roff * cosd(inpMC.angle);

% choose representative z for superficial/deep bands (midpoints from cfg)
z_sup = mean(cfg.layers.superficial_band_m) + z0;
z_dep = mean(cfg.layers.deep_band_m)        + z0;

% --- place MANY triplets (tile across width; counts from cfg)
[maskA,maskLV,maskRV] = deal(false(Nz,Nx));
maskA  = maskA  | place_triplets(cfg.vessel_counts.num_superficial_triplets, z_sup,  0,   0,  rA);
maskLV = maskLV | place_triplets(cfg.vessel_counts.num_superficial_triplets, z_sup, -dxA, dzA, rLV);
maskRV = maskRV | place_triplets(cfg.vessel_counts.num_superficial_triplets, z_sup, +dxA, dzA, rRV);

maskA  = maskA  | place_triplets(cfg.vessel_counts.num_deep_triplets, z_dep,  0,   0,  rA);
maskLV = maskLV | place_triplets(cfg.vessel_counts.num_deep_triplets, z_dep, -dxA, dzA, rLV);
maskRV = maskRV | place_triplets(cfg.vessel_counts.num_deep_triplets, z_dep, +dxA, dzA, rRV);

% keep vessels in dermis only
maskA (mask_epi | mask_water) = false;
maskLV(mask_epi | mask_water) = false;
maskRV(mask_epi | mask_water) = false;

% --- assign venous sO2 per-component (constant or layer distributions)
CC_lv = bwconncomp(maskLV,8);
CC_rv = bwconncomp(maskRV,8);
if strcmpi(cfg.optics.venous_model,'constant')
    s_lv = ones(CC_lv.NumObjects,1) * cfg.optics.S_vein_default;
    s_rv = ones(CC_rv.NumObjects,1) * cfg.optics.S_vein_default;
else
    s_lv = assign_comp_so2(CC_lv, z_grid, z_sup, z_dep, ...
        cfg.optics.superficial_vein_mean, cfg.optics.superficial_vein_std, ...
        cfg.optics.deep_vein_mean,        cfg.optics.deep_vein_std);
    s_rv = assign_comp_so2(CC_rv, z_grid, z_sup, z_dep, ...
        cfg.optics.superficial_vein_mean, cfg.optics.superficial_vein_std, ...
        cfg.optics.deep_vein_mean,        cfg.optics.deep_vein_std);
end
S_art = clamp01(normrnd(cfg.optics.artery_mean, cfg.optics.artery_std));

% --- sensor (from cfg)
col12 = round(cfg.sensor.lateral_coverage * Nx);
col12(1) = max(1,col12(1)); col12(2) = min(Nx,col12(2));
sensor.mask = false(Nz,Nx);
sensor.mask(1, col12(1):col12(2)) = true;

% --- run per wavelength: build μa, φ, p0, k-Wave forward
numSensors  = nnz(sensor.mask);
numTime     = numel(kgrid.t_array);
sensor_data = zeros(numSensors, numTime, L, 'single');
mu_a = zeros(Nz,Nx,L,'single');

for k = 1:L
    mu = mu_a_bg(:,:,k);

    % artery spectrum at this wavelength
    mua_art = cm_to_m(calc_mua(wavelength, S_art, B_blood, W_blood, F_blood, M_blood));
    mu(maskA) = single(mua_art(k));

    % paint each LV/RV component with its sO2
    mu = paint_components(mu, CC_lv, s_lv, wavelength, k, B_blood, W_blood);
    mu = paint_components(mu, CC_rv, s_rv, wavelength, k, B_blood, W_blood);
    mu_a(:,:,k) = mu;

    % φ(z,λ): lateral-mean from MC divided by background μa mean (depth-wise)
    Hmean_z = mean(single(H_mc(:,:,k)),2);          % Nz x 1
    mu_bg_z = mean(single(mu_a_bg(:,:,k)),2) + eps('single');
    phi_z   = Hmean_z ./ mu_bg_z;
    phi_z   = phi_z / max(phi_z(1),1);              % normalize at surface
    phi2D   = repmat(phi_z, [1 Nx]);                % [Nz x Nx]

    % p0
    source.p0 = mu .* phi2D;                        % Γ=1 (relative scale)

    % k-Wave forward
    sd = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
         'PlotSim', false,'DataCast','gpuArray-single'); %'DataCast', cfg.sim.datacast
    sensor_data(:,:,k) = sd;
end

% RF (time x sensors x wavelengths)
RF  = { permute(sensor_data,[2 1 3]) };

% quick PA backprojection (delay-and-sum)
PAI = { pa_recon_kwave(sensor_data, Fs, medium.sound_speed, dx, dz, sensor.mask) };

% optional US (placeholder off by default)
US = { zeros(Nz,Nx,L,'single') };

% parameters for Method A (all in millimetres)
parameters.PArfmmHeight    = medium.sound_speed * dt * 1e3;  % RF depth step (c*dt)
parameters.PAimgmmHeight   = dz * 1e3;                       % image pixel height (dz)
parameters.PApixelmmHeight = dz * 1e3;                       % <-- add this line
parameters.PApixelmmWidth  = dx * 1e3;                       % image pixel width (dx)
parameters.PAdepthOffset   = 0;
parameters.gel_mm          = cfg.layers.water_surface_z_m  * 1e3;
parameters.epidermis_mm    = cfg.layers.epidermis_thick_m  * 1e3;
parameters.tissue_mm       = (max(z_grid(:)) - cfg.layers.water_surface_z_m) * 1e3;


energy = ones(1,L,'double');
save(fullfile(save_dir,'MethodAData.mat'), ...
     'wavelength','Fs','PAI','RF','US','parameters','energy','-v7.3');
fprintf('Saved %s\n', fullfile(save_dir,'MethodAData.mat'));

% ---------------- helpers (no hardcoded numbers) ----------------
    function M = place_triplets(Ntrip, zc, dx_off, dz_off, r)
        if nargin<3, dx_off=0; end
        if nargin<4, dz_off=0; end
        if nargin<5, r = rA; end
        M = false(Nz,Nx);
        pad = 2*max([rA rLV rRV]);
        xs = linspace(pad, Lx_m - pad, max(Ntrip,2)); xs = xs(1:Ntrip);
        for i = 1:Ntrip
            x0 = xs(i);
            M = M | (((x_grid-(x0+dx_off)).^2 + (z_grid-(zc+dz_off)).^2) <= r^2);
        end
    end
end

% --------------- file-scope helpers ---------------
function y = cm_to_m(x_cm_inv), y = single(x_cm_inv * 100); end
function s = clamp01(s), s = max(0, min(1, s)); end
function [S,B,W,F,M] = unpack(v), S=v(1); B=v(2); W=v(3); F=v(4); M=v(5); end

function mu = paint_components(mu, CC, s_list, wl, k, B_blood, W_blood)
for i = 1:CC.NumObjects
    s = s_list(i);
    mu_comp = cm_to_m(calc_mua(wl, s, B_blood, W_blood, 0, 0));
    mu(CC.PixelIdxList{i}) = single(mu_comp(k));
end
end

function s_list = assign_comp_so2(CC, z_grid, z_sup, z_dep, m_sup, s_sup, m_dep, s_dep)
z_mid = 0.5*(z_sup + z_dep);
s_list = zeros(CC.NumObjects,1);
for i = 1:CC.NumObjects
    [iz,~] = ind2sub(size(z_grid), CC.PixelIdxList{i});
    z_here = mean(z_grid(iz));
    if z_here < z_mid
        s_list(i) = clamp01(normrnd(m_sup, s_sup));
    else
        s_list(i) = clamp01(normrnd(m_dep, s_dep));
    end
end
end

% function PAI = pa_backproject(sensor_data, Fs, c, dx, dz, sens_mask)
% [Ns, Nt, L] = size(sensor_data);
% cols = find(sens_mask(1,:)); x_sens = (cols-1)*dx;
% Nx = size(sens_mask,2); Nz = size(sens_mask,1);
% xv = (0:Nx-1)*dx; zv = (0:Nz-1)*dz; t = (0:Nt-1)/Fs;
% apo = hann(numel(cols)).';
% PAI = zeros(Nz, Nx, L, 'single');
% for k = 1:L
%     rf = squeeze(sensor_data(:,:,k));  % Ns x Nt
%     img = zeros(Nz,Nx,'single');
%     for ix = 1:Nx
%         for iz = 1:Nz
%             r = hypot(xv(ix) - x_sens, zv(iz));
%             tau = r / c;
%             s = zeros(1,numel(cols));
%             for sIdx = 1:numel(cols)
%                 s(sIdx) = interp1(t, rf(sIdx,:), tau(sIdx), 'linear', 0);
%             end
%             img(iz,ix) = sum(apo .* s);
%         end
%     end
%     PAI(:,:,k) = max(img,0);
% end
% end

function PAI = pa_recon_kwave(sensor_data, Fs, c, dx, dz, sens_mask)
% Fast 2D recon using kspaceLineRecon (frequency–wavenumber method).
% Accepts sensor_data as Ns x Nt  OR Ns x Nt x L. Returns Nz x Nx x L.
%
% Ns = number of sensors (active columns in sens_mask)
% Nt = time samples
% L  = number of wavelengths (pages)

% gather to CPU if needed
if isa(sensor_data,'gpuArray'), sensor_data = gather(sensor_data); end

if ndims(sensor_data) == 2
    sensor_data = reshape(sensor_data, size(sensor_data,1), size(sensor_data,2), 1);
end
[Ns, Nt, L] = size(sensor_data);

Nz = size(sens_mask,1);
Nx = size(sens_mask,2);

dt       = 1/Fs;
dx_sens  = dx;            % sensors are consecutive columns -> pitch = dx
dz_recon = c*dt/2;        % native depth sampling of the recon
z_recon  = (0:Nt-1) * dz_recon;
z_tgt    = (0:Nz-1) * dz;

PAI = zeros(Nz, Nx, L, 'single');

for k = 1:L
    sd = sensor_data(:,:,k);    % Ns x Nt
    sd_tNs = sd.';              % Nt x Ns (2-D transpose only)

    % Minimal, version-agnostic call (no optional args)
    p_xy = kspaceLineRecon(sd_tNs, dt, dx_sens, c);    % Ny x Nx_recon

    % Resample depth to our dz
    Ny_recon   = size(p_xy,1);
    z_recon_k  = z_recon(1:Ny_recon);
    PAIk       = single(interp1(z_recon_k, double(p_xy), z_tgt, 'linear', 0));  % Nz x Nx_recon

    % Fit horizontally to our Nx
    Nx_recon = size(PAIk,2);
    if Nx_recon >= Nx
        PAI(:,:,k) = PAIk(:,1:Nx);
    else
        PAI(:,:,k) = [PAIk, zeros(Nz, Nx-Nx_recon, 'single')];
    end
end
end
