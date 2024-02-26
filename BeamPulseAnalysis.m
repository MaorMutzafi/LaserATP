% Contains functions for calculating beam divergence, pulse duration, and plotting results

function [divergence_mrad, waistWidth, duration] = BeamPulseAnalysis(img, pixelSize_mm, focalLength_mm, wavelength_mm, signal, fs)
    % Main function for analyzing the beam profile and pulse duration
    
    % Analyze beam divergence and calculate waist width
    [img, divergence_mrad, waistWidth, radius_mm, beam_center] = calc_beam_div(img, pixelSize_mm, focalLength_mm, wavelength_mm);

    % Calculate pulse duration and repetition rate
    [duration, repetitionRate, pulse, t, ts, tf] = calc_pulse_dur(signal, fs);

    % Plot results
    plot_res(img, divergence_mrad, waistWidth, radius_mm, pixelSize_mm, beam_center, duration, repetitionRate, pulse, t, ts, tf)
end

function [img, divergence_mrad, waistWidth, radius_mm, beam_center] = calc_beam_div(img, pixelSize_mm, focalLength_mm, wavelength_mm)
    % Calculates beam divergence and waist width from a beam profile image

    % Convert image to grayscale and apply a filter
    img = rgb2gray(img);
    b = fir1(48, .1);
    ImgMax = max(conv2(b, b, img, 'valid'), [], 'all');
    threshold = ImgMax / 2;

    % Find the beam center and radius
    [rows, columns, ~] = size(img);
    [X, Y] = meshgrid(1:columns, 1:rows);
    centerX = mean(X(img > threshold));
    centerY = mean(Y(img > threshold));
    beam_center = [centerX, centerY];
    radius = sqrt(sum((X(img > threshold) - centerX).^2 + (Y(img > threshold) - centerY).^2) / nnz(img > threshold));

    % Convert radius to mm and calculate divergence and waist width
    radius_mm = radius * pixelSize_mm;
    divergence_rad = radius_mm / focalLength_mm;
    divergence_mrad = divergence_rad * 1000;
    waistWidth = (wavelength_mm / pi) / divergence_rad;
end

function [duration, repetitionRate, pulse, t_win, ts, tf] = calc_pulse_dur(signal, fs)
    % Calculates the pulse duration and repetition rate of a signal

    % Time vector and peak detection
    t = 0:1/fs:(length(signal)-1)/fs;
    [~, locs] = findpeaks(signal, 'MinPeakProminence', mean(signal));
    
    % Calculate repetition rate
    pulseIntervals = diff(t(locs));
    repetitionRate = mean(pulseIntervals);
    
    % Select a representative pulse for duration calculation
    middleIndex = round(length(locs)/2);
    if middleIndex == 0
        error('No pulses detected');
    end

    % Define window around the chosen pulse
    if middleIndex > 1 && middleIndex < length(locs)
        winStart = floor((locs(middleIndex-1) + locs(middleIndex)) / 2);
        winEnd = floor((locs(middleIndex) + locs(middleIndex+1)) / 2);
    else
        winStart = 1;
        winEnd = length(signal);
    end

    % Isolate the pulse and compute its duration
    pulse = signal(winStart:winEnd);
    t_win = t(winStart:winEnd);
    cumulative_energy = cumsum(pulse) / max(cumsum(pulse));
    total_energy = cumulative_energy(end);
    desired_energy = 0.9 * total_energy;
    min_interval = Inf;
    ind1 = 0;
    ind2 = 0;
    for i = 1:length(cumulative_energy)
        for j = i:length(cumulative_energy)
            if (cumulative_energy(j) - cumulative_energy(i)) >= desired_energy
                interval = j - i;
                if interval < min_interval
                    min_interval = interval;
                    ind1 = i;
                    ind2 = j;
                end
                break;
            end
        end
    end
    ts = t(winStart+ind1);
    tf = t(winStart+ind2);
    duration = tf - ts;
end

function plot_res(img, divergence_mrad, waistWidth, radius_mm, pixelSize_mm, beam_center, duration, repetitionRate, pulse, t, ts, tf)
    % Plots the results of the beam profile and pulse analysis

    % Plot signal and beam divergence
    figure;
    subplot(1, 2, 1); % Signal plot
    plot(t, pulse, 'b', 'LineWidth', 2);
    title('Pulse Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    hold on;
    plot([ts, ts], [min(pulse), max(pulse)], '--', 'color', [0, .7, 1], 'LineWidth', 2);
    plot([tf, tf], [min(pulse), max(pulse)], '--', 'color', [0, 1, .7], 'LineWidth', 2);
    hold off;
    legend({'Signal', '90% Energy Interval Starts', '90% Energy Interval Ends'});
    annotation('textbox', [0.15, 0.5, 0.3, 0.1], 'String', sprintf('Duration: %.2e sec\nRate: %.2e Hz', duration, 1/repetitionRate), 'FitBoxToText', 'on', 'BackgroundColor', 'white');

    subplot(1, 2, 2); % Beam profile plot
    x_crd = pixelSize_mm * ((0:size(img,1)-1) - beam_center(2));
    y_crd = pixelSize_mm * ((0:size(img,2)-1) - beam_center(1));
    imagesc(y_crd, x_crd, img);
    title(sprintf('Beam Profile (Divergence: %.2f mrad, Waist: %.2f mm)', divergence_mrad, waistWidth));
    xlabel('[mm]');
    ylabel('[mm]');
    hold on;
    viscircles([0 0], radius_mm, 'Color', 'r', 'LineWidth', 2.5);
    hold off;
    axis equal tight;
    colormap hot;
    text('BackgroundColor',[1 1 1], 'String', sprintf('Divergence: %.2f mrad\nWaist Width: %.2f mm\nFWHM: %.2f mm', divergence_mrad, waistWidth, radius_mm), 'Position', [-beam_center(1)*pixelSize_mm, (10-beam_center(2))*pixelSize_mm 0]);
end
