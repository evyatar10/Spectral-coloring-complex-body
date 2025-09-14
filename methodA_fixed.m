function methodA_fixed(data_mat)
% methodA_fixed(data_mat)
% Method A with:
% - RF band selection (safe 6-dB) and bandpass
% - Depth-wise fluence equalization vs ~800 nm
% - Envelope detection to PA image stack
% - NNLS unmix (Hb/HbO2) WITHOUT casting to double
% - RF plots use RF depth (c*dt); image/sO2 plots use image depth (dz)
%
% Usage:
%   methodA_fixed('./created_data_fixed/MethodAData.mat');

if nargin < 1 || ~isfile(data_mat)
    error('Provide MethodAData.mat produced by the generator');
end

S = load(data_mat);  % wavelength, Fs, PAI, RF, US, parameters, energy
required = {'wavelength','Fs','RF','parameters','energy'};
for i=1:numel(required)
    if ~isfield(S,required{i}), error('Missing %s', required{i}); end
end

% Keep native precision (single)
wavelength = S.wavelength(:)';       % 1 x L
Fs         = S.Fs;                   % scalar
RF         = S.RF{1};                % [Nt x Ns x L]
pars       = S.parameters;
energy     = S.energy;

% Ensure we have both RF-depth step and IMAGE-depth step
if ~isfield(pars,'PArfmmHeight') && isfield(pars,'PApixelmmHeight')
    pars.PArfmmHeight = pars.PApixelmmHeight;   % backward-compat (RF step)
end
if ~isfield(pars,'PAimgmmHeight')
    % fallback (if generator didn't provide dz explicitly)
    pars.PAimgmmHeight = pars.PArfmmHeight;
end
if ~isfield(pars,'PAdepthOffset'), pars.PAdepthOffset = 0; end
if ~isfield(pars,'PApixelmmWidth'), pars.PApixelmmWidth = 1; end

% 0) Spectral design matrix (keep single)
A = getHemoglobinMu(wavelength);
M = single([A.mudeoxy(:), A.muoxy(:)]);   % [L x 2]

% 1) Energy normalization
for k=1:numel(energy)
    RF(:,:,k) = RF(:,:,k) ./ max(energy(k), eps('single'));
end

% 2) SAFE 6-dB band selection
[FL, FH] = local_band_select_safe(RF, Fs);
local_plot_band_debug(RF, Fs, wavelength, FL, FH);

% 3) Bandpass
RFb = local_bandpass(RF, Fs, FL, FH);

% 4) Fluence equalization vs ~800 nm
[~,ref_idx] = min(abs(double(wavelength) - 800)); % index only
RFc = local_fluence_equalize(RFb, ref_idx);

% 5) Envelope -> PA image stack (RF-domain envelope)
PA_before = local_rftopa(RFb);   % [Nt x Ns x L]
PA_after  = local_rftopa(RFc);

% Optional overview (uses current pars scaling)
local_plot_pa_us_overview(PA_before, pars, S.US, 'PA/US overview (no ROI)', 50);

% 6) NNLS unmix to sO2 (image domain)
[Cb, SO2b] = local_unmix(PA_before, M);
[Ca, SO2a] = local_unmix(PA_after,  M);

% ---- Axes for IMAGE plots (use dz) ----
dz_mm  = pars.PAimgmmHeight;
dx_mm  = pars.PApixelmmWidth;
DepthAxis_mm = (0:size(PA_before,1)-1) * dz_mm + pars.PAdepthOffset;
WidthAxis_mm = (0:size(PA_before,2)-1) * dx_mm;

% (Optional) crop to dermis if provided
if all(isfield(pars, {'gel_mm','epidermis_mm','tissue_mm'}))
    z0 = pars.gel_mm + pars.epidermis_mm;   % start of vascular dermis
    z1 = z0 + pars.tissue_mm;               % bottom of tissue
    i0 = max(1, round(z0/dz_mm));
    i1 = min(size(SO2b,1), round(z1/dz_mm));
    SO2b = SO2b(i0:i1,:); SO2a = SO2a(i0:i1,:);
    DepthAxis_mm = DepthAxis_mm(i0:i1);
end

% ---- Plots: sO2 maps and depth curves (IMAGE depth) ----
figure('Name','sO2 maps (Method A)');
subplot(1,2,1);
imagesc(WidthAxis_mm,DepthAxis_mm,SO2b,[0 1]); axis image; set(gca,'YDir','reverse');
colorbar; title('sO2 BEFORE'); xlabel('mm'); ylabel('mm');

subplot(1,2,2);
imagesc(WidthAxis_mm,DepthAxis_mm,SO2a,[0 1]); axis image; set(gca,'YDir','reverse');
colorbar; title('sO2 AFTER');  xlabel('mm'); ylabel('mm');

figure('Name','sO2 vs depth (Method A)');
plot(DepthAxis_mm, 100*median(SO2b,2,'omitnan'), 'LineWidth',2); hold on;
plot(DepthAxis_mm, 100*median(SO2a,2,'omitnan'), 'LineWidth',2);
legend('Before','After','Location','southeast');
xlabel('Depth (mm)'); ylabel('sO2 (%)'); grid on; ylim([0 100]);

% Save results
outdir = fileparts(data_mat);
save(fullfile(outdir,'MethodA_results_fixed.mat'), 'SO2b','SO2a','Cb','Ca','wavelength','FL','FH','-v7.3');
end

% ---------- helpers (single-precision) ----------
function [FL, FH] = local_band_select_safe(RF, Fs)
% Average spectrum along depth; pick 6-dB band with guard rails.
[Nt, Ns, L] = size(RF);
X = single(0);
for k = 1:L
    R = RF(:,:,k);
    R = R - mean(R,1,'omitnan');      % per A-line DC removal
    X = X + abs(fft(R, [], 1));       % Nt x Ns
end
mag = mean(X,2,'omitnan') / max(1,(Ns*L)); mag = mag(:);
f   = (0:Nt-1).' * (Fs/Nt);
f1  = f(1:floor(Nt/2)); 
m1  = mag(1:floor(Nt/2));

th    = max(m1) * 10^(-6/20);         % 6 dB window
above = m1 >= th;
if ~any(above)
    FL = 0.5e6; 
    FH = min(12e6, Fs/2);
else
    [~,pk] = max(m1);
    l = pk; r = pk;
    while l>1 && above(l-1), l=l-1; end
    while r<numel(above) && above(r+1), r=r+1; end
    FL = f1(l); FH = f1(r);
end

% --- SAFETY RAILS ---
% FL = max(FL, 0.5e6);
% FH = min(FH, Fs/2);
% if FH <= FL*1.05
%     FH = min(max(FL*2, 2e6), Fs/2);
% end
end

function RFb = local_bandpass(RF, Fs, FL, FH)
Nt = size(RF,1);
f = linspace(-Fs/2, Fs/2, Nt).';
mask = (abs(f)>=FL & abs(f)<=FH);
mask = reshape(mask, [Nt 1 1]);
RFb = zeros(size(RF), 'like', RF);
for k=1:size(RF,3)
    X = fftshift(fft(RF(:,:,k), [], 1), 1);
    X = X .* mask;
    RFb(:,:,k) = real(ifft(ifftshift(X,1), [], 1));
end
end

function RFc = local_fluence_equalize(RFb, ref_idx)
E = squeeze(mean(RFb.^2, 2));      % Nt x L
E = max(E, eps('single'));
w = 11; if mod(w,2)==0, w=w+1; end
h = ones(w,1,'single')/w;
for k=1:size(E,2), E(:,k) = conv(E(:,k), h, 'same'); end
Eref  = E(:,ref_idx);
alpha = Eref ./ E; alpha = max(min(alpha, single(10)), single(0.1));
RFc = RFb;
for k=1:size(RFc,3)
    RFc(:,:,k) = RFb(:,:,k) .* sqrt(alpha(:,k));
end
end

function PA = local_rftopa(RF)
Nt = size(RF,1);
H = zeros(Nt,1,'single');
if mod(Nt,2)==0
    H(1)=1; H(Nt/2+1)=1; H(2:Nt/2)=2;
else
    H(1)=1; H(2:(Nt+1)/2)=2;
end
PA = zeros(size(RF), 'like', RF);
for k=1:size(PA,3)
    Z = ifft( fft(RF(:,:,k),[],1) .* H, [], 1 );
    PA(:,:,k) = abs(Z);
end
% gentle spatial smoothing
h2 = ones(3,3)/9;                       % <-- make the kernel DOUBLE
for k = 1:size(PA,3)
    PA(:,:,k) = imfilter(PA(:,:,k), h2, 'replicate');
end
end

function [C, SO2] = local_unmix(PA, M)
[m,n,L] = size(PA); K = size(M,2);
X = reshape(PA, [], L).'; 
Ccoeff = zeros(K, m*n, 'single');
for p = 1:(m*n)
    y = X(:,p);
    if all(y==0,'all')
        Ccoeff(:,p) = 0;
    else
        Ccoeff(:,p) = nnls_solve(M, y);
    end
end
C   = reshape(Ccoeff.', [m,n,K]);
Hb  = C(:,:,1); 
HbO2= C(:,:,2);
den = Hb + HbO2; den(den==0) = eps('single');
SO2 = HbO2 ./ den;
end

function x = nnls_solve(A,b)
% Single-precision projected gradient NNLS
A = single(A); b = single(b);
x = zeros(size(A,2),1,'single');
AtA = A.'*A; Atb = A.'*b;
alpha = single(1e-2); tol = single(1e-5);
for it=1:400
    g = AtA*x - Atb;
    x = max(x - alpha*g, 0);
    if norm(g,2) < tol, break; end
end
end
