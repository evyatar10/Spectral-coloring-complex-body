function local_plot_band_debug(RF, Fs, wavelength, FL, FH)
% Show power spectra for two wavelengths + their ratio (like original fig)
Nt = size(RF,1); Ns = size(RF,2); %#ok<NASGU>
L  = 2^nextpow2(Nt);
f  = (0:L-1) * (Fs/L);

% choose two wavelengths to display (match the old x1=2, x2=floor(w/2))
w = size(RF,3);
idx1 = min(2, w);
idx2 = max(1, round(w/2));

% periodogram averaged over sensors (per wavelength)
P = zeros(L, 2);
for j = 1:2
    k = (j==1)*idx1 + (j==2)*idx2;
    X = fft(RF(:,:,k), L, 1);         % FFT along time
    S = mean( abs(X).^2, 2 );         % average across sensors
    P(:,j) = S;
end

limFreq = min(30e6, Fs/2);            % keep plot within Nyquist
Pdb1 = 10*log10(P(:,1)+eps);
Pdb2 = 10*log10(P(:,2)+eps);

figure('Name','Power spectra + 6 dB band');
subplot(2,1,1);
plot(f, Pdb1, 'b'); hold on; plot(f, Pdb2, 'r');
xline(FL,'--k'); xline(FH,'--k');
xlim([0 limFreq]);
xlabel('f(Hz)'); ylabel('Power spectra (dB)');
title(sprintf('PowerSpectra and 6 dB bandwidth at wavelength : %d and %d', idx1, idx2));
grid on; hold off;

subplot(2,1,2);
plot(f, 10*log10( P(:,1)./(P(:,2)+eps) )); hold on;
xline(FL,'--k'); xline(FH,'--k');
xlim([0 limFreq]);
xlabel('f(Hz)'); ylabel(sprintf('P(%d) / P(%d) (dB)', idx1, idx2));
grid on; hold off;
end
