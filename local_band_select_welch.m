function [FL, FH, info] = local_band_select_welch(RF, Fs, opts)
% Welch-PSD based band selection for PA RF data (single precision).
% RF: [Nt x Ns x L], Fs: Hz
% opts.knee_db  : dB below peak to define the knee (default 20 dB)
% opts.min_bwHz : minimum bandwidth to enforce (default 1e6)
% opts.showPlot : true/false (default false)

if nargin < 3, opts = struct; end
if ~isfield(opts,'knee_db'),   opts.knee_db   = 20;   end
if ~isfield(opts,'min_bwHz'),  opts.min_bwHz  = 1e6;  end
if ~isfield(opts,'showPlot'),  opts.showPlot  = false;end

[Nt, Ns, L] = size(RF);

% ----- build a single long 1-D sequence by averaging across sensors & lambdas
% (removes per A-line DC; Hann window to tame leakage)
acc = zeros(Nt,1,'single');
for k = 1:L
    R = RF(:,:,k);
    R = R - mean(R,1,'omitnan');                 % DC per A-line
    w = hann(Nt,'periodic'); w = single(w(:));   % window along depth
    Rw = R .* w;                                 % (Nt x Ns)
    acc = acc + mean(Rw, 2, 'omitnan');         % average over sensors
end
acc = acc / max(1, Ns*L);                        % scale

% ----- Welch PSD (onesided)
segN = 1024; 
segN = min(segN, 2^floor(log2(Nt)));             % fit in Nt
if segN < 128, segN = min(Nt, 128); end
nover = floor(segN/2);
nfft  = max(1024, 2^nextpow2(segN));
[Pxx,f] = pwelch(acc, hann(segN,'periodic'), nover, nfft, Fs, 'onesided');
Pxx = single(Pxx(:)); f = single(f(:));

% normalize (peak = 0 dB), smooth slightly
Pxx = Pxx ./ max(Pxx + eps('single'));
Pdb = 10*log10(Pxx + eps('single'));
smN = 9; if mod(smN,2)==0, smN = smN+1; end
Pdb_s = movmean(Pdb, smN);

% ----- knee-based FH: first freq where PSD <= peak - knee_db and stays low
knee = -abs(single(opts.knee_db));
below = (Pdb_s <= knee);
FH = f(find(below,1,'first'));

% fallback: cumulative energy 95% if knee not found
if isempty(FH)
    cume = cumsum(Pxx) / max(sum(Pxx), eps('single'));
    FH = f(find(cume >= 0.95, 1, 'first'));
end

% ----- FL: low-frequency knee (first rise above knee level)
above = (Pdb_s >= knee);
FL = f(find(above,1,'first'));

% fallback: cumulative energy 5%
if isempty(FL)
    cume = cumsum(Pxx) / max(sum(Pxx), eps('single'));
    FL = f(find(cume >= 0.05, 1, 'first'));
end

% ----- guards: Nyquist cap + minimum bandwidth
nyq = single(Fs/2);
FL = max(single(0), min(FL, nyq*0.98));
FH = max(single(0), min(FH, nyq));
if FH < FL + single(opts.min_bwHz)
    FH = min(FL + single(opts.min_bwHz), nyq);
end

% Optional plot
if opts.showPlot
    figure('Name','Welch PSD band selection');
    subplot(2,1,1);
    plot(double(f), double(Pdb), 'Color',[.6 .6 1]); hold on;
    plot(double(f), double(Pdb_s), 'b'); yline(double(knee),'k--');
    xline(double(FL),'k--'); xline(double(FH),'k--');
    xlabel('f (Hz)'); ylabel('Power (dB)'); grid on;
    title('Welch PSD (blue), smoothed (dark blue)');

    subplot(2,1,2);
    cume = cumsum(Pxx)/max(sum(Pxx),eps('single'));
    plot(double(f), 10*log10(double(cume)+eps), 'b'); grid on;
    xlabel('f (Hz)'); ylabel('Cumulative energy (dB)');
    xline(double(FL),'k--'); xline(double(FH),'k--');
end

if nargout>2
    info.f = f; info.Pdb = Pdb; info.Pdb_s = Pdb_s; info.FL = FL; info.FH = FH;
end
end
