function wavelet_scalogram (u, dx)
% Plots Morlet wavelet scalogram of vector u with spacing dx

[wt, f, coi] = cwt (u, 'amor', 'VoicesPerOctave', 24, dx);    % Morlet wavelet transform in terms of frequency f

pcolor(x/1e3, 1./f/1e3, abs(wt).^2); shading flat; hold on;   % plot scalogram
plot(t, 1./coi, 'w--', 'linewidth', 2);                       % plot cone of influence for edges
set(gca,'YScale','log')
xlabel('x (km)'); ylabel('scale (km)');

% Windowed Fourier energy spectrum
[pxx, freq] = periodogram(u, hamming(length(u)), [], dx);

% Wavelet energy spectrum (wavenumber k = 2 pi f)
wav_spec = sum(abs(wt).^2, 2) * dx./f;
