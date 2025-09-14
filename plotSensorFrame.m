function plotSensorFrame(sensor_data_k, wavelength_nm, dt, c)
%PLOTSENSORFRAME  Plot k-Wave sensor RF data vs depth (mm) with zero as white.
%
%   plotSensorFrame(sensor_data_k, wavelength_nm, dt, c)
%     sensor_data_k : [Nsensors x Nt] RF array for one wavelength
%     wavelength_nm : scalar wavelength (nm) for the title
%     dt            : time step (s)
%     c             : sound speed (m/s)  (PA = one-way depth)
%
%   X-axis is depth in mm: depth = c * dt * (0:Nt-1).

    % Ensure numeric single (no need to upcast to double)
    I = single(sensor_data_k);

    % Build RF depth axis (mm), one-way
    Nt = size(I, 2);
    depth_axis_mm = (0:Nt-1) * dt * c * 1e3;

    % Robust symmetric color limits around zero
    L = prctile(abs(I(:)), 99);
    if ~isfinite(L) || L <= 0, L = max(abs(I(:))); end
    if ~isfinite(L) || L <= 0, L = 1; end

    % Plot
    figure(10); clf;
    imagesc(depth_axis_mm, 1:size(I,1), I);
    axis xy;                          % sensors top-to-bottom, depth increasing right
    xlabel('Depth (mm)');
    ylabel('Sensor index');
    title(sprintf('Sensor RF @ %g nm', wavelength_nm));
    caxis([-L L]);
    colormap(cmapZeroWhite(256));
    colorbar;
    drawnow;
end

function cmap = cmapZeroWhite(n)
%CMAPZEROWHITE  Diverging colormap with white at zero (blue↔white↔red).
    if nargin < 1, n = 256; end
    m = floor(n/2);
    neg = [linspace(0,1,m)', linspace(0,1,m)', ones(m,1)];  % Blue→White
    pos = [ones(m,1), linspace(1,0,m)', linspace(1,0,m)'];  % White→Red
    cmap = [neg; pos];
end
