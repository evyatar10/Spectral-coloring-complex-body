function [mu_a, phi] = local_buildMuAandFluence_network(wv, X, Z, mask_epi, mask_dermis, mask_vess, so2_map, gel)
% Build total mu_a (mm^-1) per wavelength and a single-source, laterally uniform fluence phi.
% - Epidermis: melanin only (no Hb)
% - Dermis/Hypodermis background: small blood volume fraction with sO2_bg
% - Vessels: whole-blood Hb with per-vessel sO2 (replaces background Hb)
% - Fluence: phi(lambda, z) = exp( - mu_eff_bg(lambda) * depth_mm ), uniform in x

[Nz, Nx] = size(X); W = numel(wv);
mu_a = zeros(Nz, Nx, W, 'single');
phi  = zeros(Nz, Nx, W, 'single');

A = getHemoglobinMu(wv);         % Hb/HbO2 (Prahl OMLC)
muHbO2 = A.muoxy(:).';           % 1 x W, mm^-1 for whole blood (150 g/L)
muHb   = A.mudeoxy(:).';

% Melanin / water / fat references (Jacques)
muMel  = A.mumel(:).';           % 1 x W
muWater= A.muwater(:).';
muFat  = A.mufat(:).';

% --- Composition parameters (simple, literature-guided) ---
% background blood volume fraction (capillary bed surrogate)
f_bv_bg   = 0.02;              % 2% blood vol. fraction in dermis
sO2_bg    = 0.8;              %tissue oxigination%%%%%%%%%%%%%%

% vessel blood volume fraction ~ whole blood
f_bv_vess = 1.0;

% bulk water/fat weightings (coarse)
f_water = 0.8;  f_fat = 0.20;

% epidermal melanin volume fraction (avascular epidermis)
f_mel_epi = 0.02;  % 2% (adjustable pigment level)

% Reduced scattering μs' per Jacques power law (skin mean)
% μs'(λ) = a * (λ/500 nm)^(-b);  a≈46 cm^-1 (4.6 mm^-1), b≈1.42  (skin mean)
a_mm = 4.6; b = 1.42;                                         % mm^-1
muspp = a_mm * (wv/500) .^ (-b);                              % 1 x W

% --- Build background absorption maps (no vessels) ---
% Epidermis (melanin only)
mu_epi = f_mel_epi * muMel + f_water*muWater + f_fat*muFat;   % 1 x W

% Dermis/Hypodermis background (Hb + water + fat)
mu_bgHb = f_bv_bg * (sO2_bg*muHbO2 + (1 - sO2_bg)*muHb);      % 1 x W
mu_bg   = mu_bgHb + f_water*muWater + f_fat*muFat;            % 1 x W

% Pre-allocate per-wavelength 2D maps
for k = 1:W
    mu_map = zeros(Nz, Nx, 'single');

    % epidermis avascular (melanin dominates)
    mu_map(mask_epi)   = single(mu_epi(k));

    % dermis background
    mu_map(mask_dermis)= single(mu_bg(k));

    % overwrite Hb term inside vessels with whole-blood Hb at per-vessel sO2
    if any(mask_vess(:))
        % Build vessel Hb contribution (replaces bg Hb)
        mu_v = f_bv_vess * ( so2_map.*muHbO2(k) + (1 - so2_map).*muHb(k) );
        % Where NaN (non-vessel), keep as 0
        mu_v(isnan(mu_v)) = 0;

        % Replace only the Hb portion: subtract bg Hb and add vessel Hb
        mu_map(mask_dermis) = mu_map(mask_dermis) - single(mu_bgHb(k));
        mu_map = mu_map + single(mu_v);
    end

    mu_a(:,:,k) = mu_map;
end

% --- Single, laterally-uniform surface illumination (no dual beams) ---
% Use background-only μa to compute μeff(λ) and φ(z) = exp(-μeff * depth_mm)
depth_mm = max((Z - gel)*1e3, 0);      % mm from tissue surface
for k = 1:W
    mu_a_bg_k = zeros(Nz, Nx, 'single');
    mu_a_bg_k(mask_epi)    = single(mu_epi(k));
    mu_a_bg_k(mask_dermis) = single(mu_bg(k));

    mu_eff_k = sqrt(3 * mu_a_bg_k .* (mu_a_bg_k + muspp(k)));  % 1/mm
    phi_k = exp( - mu_eff_k .* depth_mm );
    phi_k(~(mask_epi | mask_dermis)) = 0;   % outside tissue
    phi(:,:,k) = single(phi_k);
end
end
