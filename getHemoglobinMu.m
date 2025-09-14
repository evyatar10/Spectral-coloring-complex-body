function absorption = getHemoglobinMu(wavelength)
% getHemoglobinMu  Hemoglobin & background absorption [mm^-1] at given lambda.
% Uses Prahl's tabulated molar extinction coefficients (Excel file).
% Expected file: 'Prahl_Hb_extinction.xlsx' with columns:
%   [lambda_nm, e_HbO2_cm1_M, e_Hb_cm1_M]

fname = 'Prahl_Hb_extinction.xlsx';
if ~isfile(fname)
    error('File %s not found. Place it next to this script or update the path.', fname);
end

S = readmatrix(fname);  % columns: lambda [nm], e_HbO2 [cm^-1/M], e_Hb [cm^-1/M]
lambda_p  = S(:,1);
e_HbO2_p  = S(:,2);
e_Hb_p    = S(:,3);

% interpolate to requested wavelengths
e_HbO2 = interp1(lambda_p, e_HbO2_p, wavelength, 'pchip', 'extrap');
e_Hb   = interp1(lambda_p, e_Hb_p,   wavelength, 'pchip', 'extrap');

% convert extinction -> absorption mu [cm^-1]
C_Hb = 150;      % g/L (whole-blood Hb concentration)
MW   = 64500;    % g/mol
mu_cm_HbO2 = 2.303 .* e_HbO2 .* C_Hb ./ MW;
mu_cm_Hb   = 2.303 .* e_Hb   .* C_Hb ./ MW;

% convert cm^-1 to mm^-1
mu_mm_HbO2 = mu_cm_HbO2 / 10;
mu_mm_Hb   = mu_cm_Hb   / 10;

absorption.muoxy   = mu_mm_HbO2(:).';
absorption.mudeoxy = mu_mm_Hb(:).';

  %% 2) Melanin (Jacques model)
  % µa(cm^-1) = 1.70e12·λ^-3.48  → mm^-1
  absorption.mumel = 1.70e12 .* wavelength.^(-3.48) * 0.1;

  %% 3) Fat (van Veen et al. JBO 10(5) 054004)
  lam_fat    = [760, 930];      % nm
  mu_fat_m1  = [1.73, 12.8];    % m^-1
  mu_fat_mm1 = mu_fat_m1 * 0.001;  % mm^-1

  % interpolate within range
  absorption.mufat = interp1(lam_fat, mu_fat_mm1, wavelength, 'pchip');
  % clamp to endpoint values outside data range
  absorption.mufat(wavelength < lam_fat(1)) = mu_fat_mm1(1);
  absorption.mufat(wavelength > lam_fat(end)) = mu_fat_mm1(end);

  %% 4) Water (Hale & Querry Appl. Opt. 12(3) 555–563)
  lam_w    = [825,850,875,900,925,950,975,1000];  % nm
  mu_w_mm1 = [0.0282,0.0433,0.0562,0.0679,0.144,0.388,0.449,0.36];  % mm^-1

  absorption.muwater = interp1(lam_w, mu_w_mm1, wavelength, 'pchip');
  % clamp to endpoint values outside data range
  absorption.muwater(wavelength < lam_w(1))   = mu_w_mm1(1);
  absorption.muwater(wavelength > lam_w(end)) = mu_w_mm1(end);

% Minimal background terms (kept for shape compatibility)
% absorption.mumel   = zeros(size(wavelength));
% absorption.mufat   = linspace(0.002, 0.014, numel(wavelength));
% absorption.muwater = max(0, interp1([825 950 975 1000],[0.028 0.388 0.449 0.36], wavelength, 'pchip','extrap'));
end
