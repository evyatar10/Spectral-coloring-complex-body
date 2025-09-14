function local_plot_pa_us_overview(PA_cube, pars, us_in, fig_name, dr_db)
% local_plot_pa_us_overview  Show PA (dB) + US overview (no ROI).
%
% PA_cube : [Nt x Ns x L] PA amplitude stack (e.g., PA_before or PA_after)
% pars    : struct with fields PApixelmmHeight, PApixelmmWidth, PAdepthOffset
% us_in   : []  OR  cell {UScube} (Nt x Ns x L)  OR  US image (Nt x Ns)
% fig_name: (optional) figure title
% dr_db   : (optional) PA dynamic range in dB (default 50)
%
% Example call (inside methodA_fixed after envelope step):
%   local_plot_pa_us_overview(PA_before, pars, S.US, 'PA/US overview (no ROI)', 50);

if nargin < 4 || isempty(fig_name), fig_name = 'PA/US overview (no ROI)'; end
if nargin < 5 || isempty(dr_db),   dr_db   = 50; end

% Axes (mm)
DepthAxis = (0:size(PA_cube,1)-1) * pars.PApixelmmHeight + pars.PAdepthOffset;
WidthAxis = (0:size(PA_cube,2)-1) * pars.PApixelmmWidth;

% PA image: mean across wavelengths, dB display
PAimg = mean(PA_cube, 3);
PAimg = max(PAimg, eps);                             % avoid log(0)
PAimg_db = 20*log10( PAimg / max(PAimg(:)) );

% US image handling (accept [] / cell / matrix)
USimg = [];
if ~isempty(us_in)
    if iscell(us_in) && ~isempty(us_in{1})
        U = us_in{1};
        if ndims(U) == 3, USimg = mean(U,3);
        elseif ismatrix(U), USimg = U;
        end
    elseif ismatrix(us_in)
        USimg = us_in;
    end
end

figure('Name', fig_name);
% PA
subplot(2,1,1);
imagesc(WidthAxis, DepthAxis, PAimg_db, [-dr_db 0]); 
axis image; set(gca,'YDir','reverse');
colormap(gca, hot); colorbar; 
title('Photoacoustic image'); xlabel('(mm)'); ylabel('(mm)');

% US (or fallback to PA if not provided)
subplot(2,1,2);
if ~isempty(USimg)
    imagesc(WidthAxis, DepthAxis, USimg);
    colormap(gca, gray);
else
    imagesc(WidthAxis, DepthAxis, PAimg_db, [-dr_db 0]);
    colormap(gca, gray);
end
axis image; set(gca,'YDir','reverse'); colorbar;
title('Ultrasound image'); xlabel('(mm)'); ylabel('(mm)');
end
