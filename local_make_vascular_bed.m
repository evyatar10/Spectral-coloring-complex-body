function [mask_vess, so2_map] = local_make_vascular_bed(X, Z, gel, opts)
% Build a 2-plexus vascular bed with many vessels, laterally staggered.
% Returns:
%   mask_vess : logical mask of all vessels
%   so2_map   : per-pixel sO2 (vessel pixels get artery/vein values; elsewhere NaN)

[Nz, Nx] = size(X);
mask_vess = false(Nz, Nx);
so2_map   = nan(Nz, Nx);

rng(0); % reproducible

% Superficial plexus (small radii)
N1 = opts.num_superficial;
r1 = linspace(opts.superficial_r_mm(1), opts.superficial_r_mm(2), N1) * 1e-3;
z1 = gel + (opts.superficial_z_mm(1) + (opts.superficial_z_mm(2)-opts.superficial_z_mm(1))*rand(1,N1)) * 1e-3;
x1 = linspace(min(X(:)), max(X(:)), N1) + (0.2e-3)*randn(1,N1); % stagger laterally

% Deep plexus (larger vessels)
N2 = opts.num_deep;
r2 = linspace(opts.deep_r_mm(1), opts.deep_r_mm(2), N2) * 1e-3;
z2 = gel + (opts.deep_z_mm(1) + (opts.deep_z_mm(2)-opts.deep_z_mm(1))*rand(1,N2)) * 1e-3;
x2 = linspace(min(X(:))+0.5e-3, max(X(:))-0.5e-3, N2) + (0.3e-3)*randn(1,N2);

% Combine centers/radii
xc = [x1 x2]; zc = [z1 z2]; rc = [r1 r2];
N  = numel(rc);

% Assign artery/vein labels and sO2
is_artery = false(1,N);
is_artery( randperm(N, round(opts.artery_fraction*N)) ) = true;
sO2v = zeros(1,N);
sO2v( is_artery) = opts.artery_so2;
sO2v(~is_artery) = opts.vein_so2;

% Rasterize circles
for i = 1:N
    mask_i = ((X - xc(i)).^2 + (Z - zc(i)).^2) <= rc(i)^2;
    mask_vess = mask_vess | mask_i;
    so2_map(mask_i) = sO2v(i);
end
end
