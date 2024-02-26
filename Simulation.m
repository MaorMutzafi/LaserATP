
% --- Simulation Script ---
% Load the beam profile image
img = imread('beam_profile.jpg');

% Define pixel size, lens properties, and laser properties
pixelSize_mm = .01;
focalLength_mm = 1000;
wavelength_mm = 1.55;

% Define sampling frequency and create a time vector
fs = 100000;
t = 0:1/fs:.1-1/fs;

% Define pulse properties and create a pulse signal
pulse_width = 0.00005;
pulse_rate = 1e3;
pulse_amplitude = 1;
signal = zeros(size(t));
t_p = rand(1)/pulse_rate;
while t_p < t(end)
    signal = signal + pulse_amplitude * exp(-(t - t_p).^2 / pulse_width^2 / 2);
    t_p = t_p + 1 / pulse_rate;
end

% Perform analysis
[divergence_mrad, waistWidth, duration] = BeamPulseAnalysis(img, pixelSize_mm, focalLength_mm, wavelength_mm, signal, fs);
