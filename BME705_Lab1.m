%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME705: Rehabilitation Engineering
% Lab 1: Applications of FES and EMG in Rehabilitation
%
% Created by: Devon Santillo, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Names: Tanvir
 
% Student IDs: 501104056
%

close all
clear all

%% Initializations
% Load all relevent files and, if necessary, save the components of them
% into variables

Data = load("TA_data1.mat");
fatigue_data = load("FESfatigue_data10.mat");
stimulation_data = load("FStim_data.mat");

%define the three sampling frequencies for each data set
fs_Data = 2000;
fs_fatigue_data = 1000;
fs_simulation_data = 100;

%% Part 1: Introduction to EMG analysis
% separate EMG and force data using the dot operator
EMG = (Data.increase_ta_emg);
Force = (Data.increase_ta_force);

% Define time as a vector

t = (0:length(EMG)-1) / fs_Data;

% Processing:

G = 500;

EMG= EMG/G;
EMG= EMG*1000; % convert to mV

%1)original graphs

figure;
plot(t, EMG);
title('EMG');
xlabel('Time (s)');
ylabel('EMG (mV)');
grid;

% 2) normalizations

%3) rectification; 

EMG_rectified = abs(EMG);
figure;
plot(t, EMG_rectified);
title('Rectified EMG');
xlabel('Time (s)');
ylabel('Rectified EMG (mV)');
grid;

%4) Filtered EMG signal
fc = 2.5; % corner frequency
%establishing the transfer function of a 4th order butterworth filter

[b, a] = butter(4, fc/(fs_Data/2));
EMG_filtered = filter(b, a, EMG_rectified);

figure;
plot(t, EMG_filtered);
title('Filtered EMG')
xlabel('Time (s)');
ylabel('EMG (mV)');
grid;

%5) calculating iEMG
iEMG = cumtrapz(t, EMG_rectified);

figure;
plot(t, iEMG);
title('iEMG')
xlabel('Time (s)');
ylabel('iEMG (mV)');
grid;

%6) Dividing signal into individual contractions

EMGr = reshape(EMG_filtered, [30000, 10]);
Forcer = reshape(Force, [30000, 10]);
tr = reshape(t, [30000, 10]);

RMS = zeros(10, 1);
average_force = zeros(10, 1);
contraction_indices = zeros(10, 2);

for i = 1:10
    limit = 0.2*max(Forcer(:, i));  % Define limit for force
    contractions = find(Forcer(:, i) > limit);  % Find contractions based on force signal

    contraction_indices(i, 1) = contractions(1);
    contraction_indices(i, 2) = contractions(end);
    RMS(i) = rms(EMGr(contraction_indices(i, 1):contraction_indices(i, 2), i));
    average_force(i) = mean(Forcer(contraction_indices(i, 1):contraction_indices(i, 2), i));
    
    % Plot the EMG signal and force for the current segment
    figure;
    
    % Plot EMG signal
    subplot(2, 1, 1);
    plot(tr(:, i), EMGr(:, i), 'b'); % Plot the entire EMG signal in blue
    hold on;
    plot(tr(contraction_indices(i, 1):contraction_indices(i, 2), i), ...
        EMGr(contraction_indices(i, 1):contraction_indices(i, 2), i), 'r', 'LineWidth', 1.5); % Plot the contraction region in red
    yline(RMS(i), '--g', 'LineWidth', 1.5); % Add a horizontal line to indicate the RMS value
    xlabel('Time (s)');
    ylabel('EMG Signal');
    title(['Segment ', num2str(i), ' - EMG Signal with Contraction and RMS Value']);
    legend('EMG Signal', 'Contraction Region', 'RMS Value');
    grid on;
    hold off;

    % Plot force signal
    subplot(2, 1, 2);
    plot(tr(:, i), Forcer(:, i), 'b'); % Plot the entire force signal in blue
    hold on;
    plot(tr(contraction_indices(i, 1):contraction_indices(i, 2), i), ...
        Forcer(contraction_indices(i, 1):contraction_indices(i, 2), i), 'r', 'LineWidth', 1.5); % Plot the contraction region in red
    yline(average_force(i), '--g', 'LineWidth', 1.5); % Add a horizontal line to indicate the average force
    xlabel('Time (s)');
    ylabel('Force');
    title(['Segment ', num2str(i), ' - Force Signal with Contraction and Average Force']);
    legend('Force Signal', 'Contraction Region', 'Average Force');
    grid on;
    hold off;
end

figure;
scatter(average_force, RMS, 'filled');
xlabel('Average Force (N)');
ylabel('RMS of EMG Signal');
title('Relationship between Average Force and RMS of EMG Signal');
grid on;

%% Part 2: Effects of stimulation frequency on muscle fatigue

% Find the shortest length for 30Hz data
min_length_30 = min([length(fatigue_data.force30_1), length(fatigue_data.force30_2), length(fatigue_data.force30_3)]);

% Trim all 30Hz trials to the shortest length
force30_1_trimmed = fatigue_data.force30_1(1:min_length_30);
force30_2_trimmed = fatigue_data.force30_2(1:min_length_30);
force30_3_trimmed = fatigue_data.force30_3(1:min_length_30);

tf30 = (0 : min_length_30 - 1) / fs_fatigue_data;

% Calculate the average and standard deviation for 30Hz trials
avg_force30 = (force30_1_trimmed + force30_2_trimmed + force30_3_trimmed) / 3;

% Plot the 30Hz data
figure;
plot(tf30, force30_1_trimmed, 'b', 'DisplayName', 'Trial 1');
hold on;
plot(tf30, force30_2_trimmed, 'r', 'DisplayName', 'Trial 2');
plot(tf30, force30_3_trimmed, 'g', 'DisplayName', 'Trial 3');
plot(tf30, avg_force30, 'k', 'LineWidth', 1.5, 'DisplayName', 'Average Force');
hold off;

xlabel('Time (s)');
ylabel('Force (N)');
title('30 Hz Stimulation Signals with Synchronized Average');
legend('show');
grid on;

% Repeat the same process for 60Hz data
min_length_60 = min([length(fatigue_data.force60_1), length(fatigue_data.force60_2), length(fatigue_data.force60_3)]);

% Trim all 60Hz trials to the shortest length
force60_1_trimmed = fatigue_data.force60_1(1:min_length_60);
force60_2_trimmed = fatigue_data.force60_2(1:min_length_60);
force60_3_trimmed = fatigue_data.force60_3(1:min_length_60);

% Create time vectors for the trimmed data
tf60 = (0 : min_length_60 - 1) / fs_fatigue_data;

% Calculate the average and standard deviation for 60Hz trials
avg_force60 = (force60_1_trimmed + force60_2_trimmed + force60_3_trimmed) / 3;

% Plot the 60Hz data
figure;
plot(tf60, force60_1_trimmed, 'b', 'DisplayName', 'Trial 1');
hold on;
plot(tf60, force60_2_trimmed, 'r', 'DisplayName', 'Trial 2');
plot(tf60, force60_3_trimmed, 'g', 'DisplayName', 'Trial 3');
plot(tf60, avg_force60, 'k', 'LineWidth', 1.5, 'DisplayName', 'Average Force');
hold off;

xlabel('Time (s)');
ylabel('Force (N)');
title('60 Hz Stimulation Signals with Synchronized Average');
legend('show');
grid on;

%% Part 3: Stimulation input - motor threshold investigation

% Define time of signal
t_stimulation_data = (0 : length(stimulation_data.force) - 1) / fs_simulation_data;

% Plot original force signal
figure;
plot(t_stimulation_data, stimulation_data.force, 'DisplayName', 'Data');

xlabel('Time (s)');
ylabel('Force (N)');
title('Raw Force Data');
legend('show');
grid on;

 
% Zero-mA Force: Subtract the dc-offset at 0mA stimulation
amount_of_points = 0;
offset_total = 0;
force_offset = stimulation_data.force;
i = 1;
start_i = 1;
while i ~= length(t_stimulation_data)
    if stimulation_data.stim_train(i) == 0
        amount_of_points = amount_of_points + 1;
        offset_total = offset_total + stimulation_data.force(i);
        i = i + 1;
    else
        offset = offset_total/amount_of_points;
        while stimulation_data.stim_train(i) ~= 0
            i = i + 1;
        end
        force_offset(start_i:i-1) = force_offset(start_i:i-1) - offset;
        amount_of_points = 0;
        offset_total = 0;
        start_i = i;
    end
end


% Plot the force and stimulation curves on the same chart. 
figure;
plot(t_stimulation_data, stimulation_data.stim_train, 'DisplayName', 'Stimulation (mA)');
hold on;
plot(t_stimulation_data, force_offset, 'DisplayName', 'Force (N)', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('Force and Stimulation Curves');
legend('show');
grid on;

%% Part 3: Stimulation input - motor threshold investigation

% Define time of signal
t_stimulation_data = (0 : length(stimulation_data.force) - 1) / fs_simulation_data;

% Plot original force signal
figure;
plot(t_stimulation_data, stimulation_data.force, 'DisplayName', 'Data');

xlabel('Time (s)');
ylabel('Force (N)');
title('Raw Force Data');
legend('show');
grid on;

 
% Zero-mA Force: Subtract the dc-offset at 0mA stimulation
amount_of_points = 0;
offset_total = 0;
force_offset = stimulation_data.force;
i = 1;
while stimulation_data.stim_train(i) == 0
    amount_of_points = amount_of_points + 1;
    offset_total = offset_total + stimulation_data.force(i);
    i = i + 1;
end
offset = offset_total/amount_of_points;
force_offset = force_offset - offset;

% Plot the force and stimulation curves on the same chart. 
figure;
plot(t_stimulation_data, stimulation_data.stim_train, 'DisplayName', 'Stimulation (mA)');
hold on;
plot(t_stimulation_data, force_offset, 'DisplayName', 'Force (N)', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('Force and Stimulation Curves');
legend('show');
grid on;