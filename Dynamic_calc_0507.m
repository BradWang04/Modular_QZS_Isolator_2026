%% Static displacement - QZS analysis
clear all;clc;
%%
a = 50; %mm
b = 80; %mm
g = 9.8;
x0 = 5;
phi_0 = asin((b-x0)/b);

k = (300/1E3*g)/(80); % N/mm

% x1 = [0:1:45]; % Displacement
phi_1 = [phi_0:-0.05:asin((b-50)/b)];

EF_0 = sqrt(3)*b*cos(phi_0)*ones(size(phi_1));
EF_1 = sqrt(3)*b*cos(phi_1);

disp = 2*b*ones(size(phi_1))-2*b*sin(phi_1);
F_Ver = 9/2*k*b*sin(phi_1).*(cos(phi_1) - cos(phi_0)*ones(size(phi_1)) )./cos(phi_1);


%----- QZS region ------
x_1 = disp(end:-1:1)-2*x0;
y_1 = F_Ver(end:-1:1);

% plot(x_1,y_1)

% 使用数值微分方法对 y_1 关于 x_1 求导
dx = diff(x_1);
dy = diff(y_1);
k_2 = dy ./ dx;

% 找出满足 abs(k_2) <= 0.5*k 的 x_1 值
valid_indices = abs(k_2) <= 0.5*k;
valid_x_1 = x_1(1:end-1); % 注意 k_2 的长度比 x_1 少 1
valid_x_1 = valid_x_1(valid_indices);
valid_y_1 = y_1(1:end-1); % 注意 k_2 的长度比 y_1 少 1
valid_y_1 = valid_y_1(valid_indices);

% 计算对应的 x_1 长度
total_x_1_length = sum(diff(valid_x_1));

% 输出结果
fprintf('满足 abs(k_2) <= 0.5*k 的 x_1 长度为: %f mm\n', total_x_1_length);

% 作图
figure()
hold on
plot(x_1, y_1, 'b', 'LineWidth', 1.5)
plot(valid_x_1, valid_y_1, 'r', 'LineWidth', 1.5)

% 计算黄色矩形的上下界和宽度
y_min = min(valid_y_1);
y_max = max(valid_y_1);
x_min = min(valid_x_1);
x_max = max(valid_x_1);

% 画出黄色、50% 透明度的矩形
rectangle('Position', [x_min, y_min, x_max-x_min, y_max-y_min], 'FaceColor', [1, 1, 0, 0.5], 'EdgeColor', 'none')

grid on
title('Mathematical Model for Force-Displacement Curve', 'FontSize', 12)
xlabel('Displacement [mm]', 'FontSize', 15)
ylabel('Loading [N]', 'FontSize', 15)
hold off


%% Dynamics
% System parameters
M = max(F_Ver)/g;
kh = k;
fa_square = 5*0.5^2;

% Frequency range with finer resolution for better curve
f_values = linspace(1, 90, 180); % Doubled resolution
za_all = sqrt(fa_square./f_values);
Tr_values = zeros(size(f_values));

% Time parameters with increased simulation time and finer steps
t_end = 10;  % Longer simulation time for better steady state
dt = 0.01;  % Smaller time step
tspan_1 = 0:dt:t_end;
phi0 = [phi_0; 0];

% ODE solver options for better accuracy
options = odeset('RelTol', 1e-8, ...
                'AbsTol', 1e-10, ...
                'MaxStep', dt);

% Progress bar initialization
fprintf('Computing transmissibility...\n');
progress_bar = waitbar(0, 'Starting calculations...');

for i = 1:length(f_values)
    try
        za = za_all(i);
        f = f_values(i);
        
        % Solve ODE with improved accuracy
        [t, phi_dyn] = ode45(@(t,phi) myODE_dyn_1(t,phi,b,phi_0,g,f,za,kh,M), ...
                            tspan_1, phi0, options);
        
        % Initialize arrays
        phi_all = zeros(length(phi_dyn), 3);
        epsilon = 1e-10; % Numerical stability factor
        
        % Calculate derivatives with improved numerical stability
        phi_ddot = 1./(2*b*cos(phi_dyn(:, 1)) + epsilon) .* ...
            (2*b*sin(phi_dyn(:, 1)).*phi_dyn(:, 2).^2 + ...
             4*pi^2*f^2*za*sin(2*pi.*f.*t) - g + ...
             (9*kh*b)/(2*M)*tan(phi_dyn(:, 1)).*(cos(phi_dyn(:, 1)) - cos(phi_0)));
        
        % Store all phi components
        phi_all(:,1) = phi_dyn(:, 1);  % phi
        phi_all(:,2) = phi_dyn(:, 2);  % phi_dot
        phi_all(:,3) = phi_ddot;       % phi_ddot
        
        % Calculate acceleration with improved stability
        z_ddot = -g + 1./(2*M*b*cos(phi_all(:,1)) + epsilon) .* ...
            (3/2*kh*b^2*sin(phi_all(:,1)).*cos(phi_all(:,1)) + ...
             4*M*b^2*sin(phi_all(:,1)).*cos(phi_all(:,1)).*phi_all(:,2).^2 - ...
             4*M*b^2*cos(phi_all(:,1)).^2.*phi_all(:,3));
        
        % Steady state analysis (last 30% of simulation)
        steady_state_start = floor(0.7 * length(t));
        t_steady = t(steady_state_start:end);
        z_ddot_steady = z_ddot(steady_state_start:end);
        phi_steady = phi_all(steady_state_start:end,:);
        
        % Calculate RMS values for transmissibility
        response = z_ddot_steady + ...
                  2*b*cos(phi_steady(:,1)).*phi_steady(:,3) - ...
                  2*b*sin(phi_steady(:,1)).*phi_steady(:,2).^2;
        
        response_rms = rms(response);
        excitation_rms = rms(z_ddot_steady);
        
        % Calculate transmissibility in dB
        if excitation_rms > 0 && response_rms > 0
            Tr = response_rms / excitation_rms;
            Tr_values(i) = 20*log10(Tr);
        else
            Tr_values(i) = NaN;
        end
        
        % Update progress bar
        waitbar(i/length(f_values), progress_bar, ...
                sprintf('Processing frequency %0.1f Hz (%d/%d)', ...
                f, i, length(f_values)));
        
    catch ME
        warning('Error at frequency %d Hz: %s', f, ME.message);
        Tr_values(i) = NaN;
    end
end

% Close progress bar
close(progress_bar);

% Remove any NaN values for plotting
valid_indices = ~isnan(Tr_values);
f_valid = f_values(valid_indices);
Tr_valid = Tr_values(valid_indices);

% Plotting with improved visualization
figure('Position', [100, 100, 800, 600])
semilogx(f_valid, Tr_valid, 'b-', 'LineWidth', 2)
hold on
grid on

% Add reference lines
yline(0, '--k', 'LineWidth', 1);
yline(-3, ':k', 'LineWidth', 1, 'Label', '-3 dB');

% Enhance plot appearance
xlabel('Frequency [Hz]', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Transmissibility [dB]', 'FontSize', 12, 'FontWeight', 'bold')
title('Frequency Response: Transmissibility vs. Frequency', ...
      'FontSize', 14, 'FontWeight', 'bold')

% Add axis limits and formatting
xlim([min(f_valid), max(f_valid)])
ylim([min(Tr_valid)-5, max(Tr_valid)+5])
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(gca, 'FontSize', 10)

% Add legend
legend('Transmissibility', 'Location', 'best')

% Save results
results.frequency = f_valid;
results.transmissibility = Tr_valid;
results.parameters = struct('M', M, 'kh', kh, 'fa_square', fa_square);

% Print summary statistics
fprintf('\nAnalysis Complete:\n');
fprintf('Frequency range: %.1f to %.1f Hz\n', min(f_valid), max(f_valid));
fprintf('Maximum transmissibility: %.2f dB\n', max(Tr_valid));
fprintf('Minimum transmissibility: %.2f dB\n', min(Tr_valid));

%% Dynamics - 2: Using local maxima version

% System parameters
M = max(F_Ver)/g;
kh = k;
fa_square = 5*0.5^2;

% Frequency range with finer resolution for better curve
f_values = linspace(1, 90, 180); % Doubled resolution
za_all = sqrt(fa_square./f_values);
Tr_values = zeros(size(f_values));

% Time parameters with increased simulation time and finer steps
t_end = 10;  % Longer simulation time for better steady state
dt = 0.001;  % Smaller time step
tspan_1 = 0:dt:t_end;
phi0 = [phi_0; 0];

% ODE solver options for better accuracy
options = odeset('RelTol', 1e-8, ...
                'AbsTol', 1e-10, ...
                'MaxStep', dt);

% Progress bar initialization
fprintf('Computing transmissibility...\n');
progress_bar = waitbar(0, 'Starting calculations...');

for i = 1:length(f_values)
    try
        za = za_all(i);
        f = f_values(i);
        
        % Solve ODE with improved accuracy
        [t, phi_dyn] = ode45(@(t,phi) myODE_dyn_1(t,phi,b,phi_0,g,f,za,kh,M), ...
                            tspan_1, phi0, options);
        
        % Initialize arrays
        phi_all = zeros(length(phi_dyn), 3);
        epsilon = 1e-10; % Numerical stability factor
        
        % Calculate derivatives with improved numerical stability
        phi_ddot = 1./(2*b*cos(phi_dyn(:, 1)) + epsilon) .* ...
            (2*b*sin(phi_dyn(:, 1)).*phi_dyn(:, 2).^2 + ...
             4*pi^2*f^2*za*sin(2*pi.*f.*t) - g + ...
             (9*kh*b)/(2*M)*tan(phi_dyn(:, 1)).*(cos(phi_dyn(:, 1)) - cos(phi_0)));
        
        % Store all phi components
        phi_all(:,1) = phi_dyn(:, 1);  % phi
        phi_all(:,2) = phi_dyn(:, 2);  % phi_dot
        phi_all(:,3) = phi_ddot;       % phi_ddot
        
        % Calculate acceleration with improved stability
        z_ddot = -g + 1./(2*M*b*cos(phi_all(:,1)) + epsilon) .* ...
            (3/2*kh*b^2*sin(phi_all(:,1)).*cos(phi_all(:,1)) + ...
             4*M*b^2*sin(phi_all(:,1)).*cos(phi_all(:,1)).*phi_all(:,2).^2 - ...
             4*M*b^2*cos(phi_all(:,1)).^2.*phi_all(:,3));
        
        % Steady state analysis (last 30% of simulation)
        steady_state_start = floor(0.7 * length(t));
        t_steady = t(steady_state_start:end);
        z_ddot_steady = z_ddot(steady_state_start:end);
        phi_steady = phi_all(steady_state_start:end,:);
        
        % Calculate top plate acceleration (response)
        response = z_ddot_steady + ...
                  2*b*cos(phi_steady(:,1)).*phi_steady(:,3) - ...
                  2*b*sin(phi_steady(:,1)).*phi_steady(:,2).^2;
        
        % Calculate bottom plate acceleration (excitation)
        excitation = 4*pi^2*f^2*za*sin(2*pi*f*t_steady);
        
        % Find local maxima for both response and excitation
        [peaks_response, ~] = findpeaks(response);
        [valleys_response, ~] = findpeaks(-response);
        valleys_response = -valleys_response;
        
        [peaks_excitation, ~] = findpeaks(excitation);
        [valleys_excitation, ~] = findpeaks(-excitation);
        valleys_excitation = -valleys_excitation;
        
        % Calculate average of absolute peak values
        avg_response = mean([abs(peaks_response); abs(valleys_response)]);
        avg_excitation = mean([abs(peaks_excitation); abs(valleys_excitation)]);
        
        % Calculate transmissibility in dB
        if avg_excitation > 0 && avg_response > 0
            Tr = avg_response / avg_excitation;
            Tr_values(i) = 20*log10(Tr);
        else
            Tr_values(i) = NaN;
        end
        
        % Update progress bar
        waitbar(i/length(f_values), progress_bar, ...
                sprintf('Processing frequency %0.1f Hz (%d/%d)', ...
                f, i, length(f_values)));
        
    catch ME
        warning('Error at frequency %d Hz: %s', f, ME.message);
        Tr_values(i) = NaN;
    end
end

% Close progress bar
close(progress_bar);

% Remove any NaN values for plotting
valid_indices = ~isnan(Tr_values);
f_valid = f_values(valid_indices);
Tr_valid = Tr_values(valid_indices);

% Plotting with improved visualization
figure('Position', [100, 100, 800, 600])
semilogx(f_valid, Tr_valid, 'b-', 'LineWidth', 2)
hold on
grid on

% Add reference lines
yline(0, '--k', 'LineWidth', 1);
yline(-3, ':k', 'LineWidth', 1, 'Label', '-3 dB');

% Enhance plot appearance
xlabel('Frequency [Hz]', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Transmissibility [dB]', 'FontSize', 12, 'FontWeight', 'bold')
title('Frequency Response: Transmissibility vs. Frequency', ...
      'FontSize', 14, 'FontWeight', 'bold')

% Add axis limits and formatting
xlim([min(f_valid), max(f_valid)])
ylim([min(Tr_valid)-5, max(Tr_valid)+5])
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(gca, 'FontSize', 10)

% Add legend
legend('Transmissibility', 'Location', 'best')

% Save results
results.frequency = f_valid;
results.transmissibility = Tr_valid;
results.parameters = struct('M', M, 'kh', kh, 'fa_square', fa_square);

% Print summary statistics
fprintf('\nAnalysis Complete:\n');
fprintf('Frequency range: %.1f to %.1f Hz\n', min(f_valid), max(f_valid));
fprintf('Maximum transmissibility: %.2f dB\n', max(Tr_valid));
fprintf('Minimum transmissibility: %.2f dB\n', min(Tr_valid));

%% Dynamics: 3
% System parameters
M = max(F_Ver)/g;
kh = k;
fa_square = 5*0.5^2;

% Frequency range with finer resolution
f_values = linspace(1, 90, 180); % Higher resolution
za_all = sqrt(fa_square./f_values);
Tr_values = zeros(size(f_values));

% Time parameters with increased simulation time
t_end = 10;  % Longer simulation time
dt = 0.01;   % Smaller time step
tspan_1 = 0:dt:t_end;
phi0 = [phi_0; 0];

% ODE solver options for better accuracy
options = odeset('RelTol', 1e-8, ...
                'AbsTol', 1e-10, ...
                'MaxStep', dt);

% Progress bar initialization
fprintf('Computing transmissibility...\n');
progress_bar = waitbar(0, 'Starting calculations...');

for i = 1:length(f_values)
    try
        za = za_all(i);
        f = f_values(i);
        
        % Solve ODE with improved accuracy
        [t, phi_dyn] = ode45(@(t,phi) myODE_dyn_1(t,phi,b,phi_0,g,f,za,kh,M), ...
                            tspan_1, phi0, options);
        
        % Initialize arrays
        phi_all = zeros(length(phi_dyn), 3);
        epsilon = 1e-10; % Numerical stability factor
        
        % Calculate derivatives
        phi_ddot = 1./(2*b*cos(phi_dyn(:, 1)) + epsilon) .* ...
            (2*b*sin(phi_dyn(:, 1)).*phi_dyn(:, 2).^2 + ...
             4*pi^2*f^2*za*sin(2*pi.*f.*t) - g + ...
             (9*kh*b)/(2*M)*tan(phi_dyn(:, 1)).*(cos(phi_dyn(:, 1)) - cos(phi_0)));
        
        % Store all phi components
        phi_all(:,1) = phi_dyn(:, 1);  % phi
        phi_all(:,2) = phi_dyn(:, 2);  % phi_dot
        phi_all(:,3) = phi_ddot;       % phi_ddot
        
        % Calculate acceleration
        z_ddot = -g + 1./(2*M*b*cos(phi_all(:,1)) + epsilon) .* ...
            (3/2*kh*b^2*sin(phi_all(:,1)).*cos(phi_all(:,1)) + ...
             4*M*b^2*sin(phi_all(:,1)).*cos(phi_all(:,1)).*phi_all(:,2).^2 - ...
             4*M*b^2*cos(phi_all(:,1)).^2.*phi_all(:,3));
        
        % Steady state analysis (last 30% of simulation)
        steady_state_start = floor(0.7 * length(t));
        t_steady = t(steady_state_start:end);
        z_ddot_steady = z_ddot(steady_state_start:end);
        
        % Calculate response acceleration
        response_accel = z_ddot_steady + ...
                       2*b*cos(phi_all(steady_state_start:end,1)).*phi_all(steady_state_start:end,3) - ...
                       2*b*sin(phi_all(steady_state_start:end,1)).*phi_all(steady_state_start:end,2).^2;
        
        % Find local maxima for excitation (z_ddot_steady)
        [exc_pks_top, exc_locs_top] = findpeaks(z_ddot_steady);
        [exc_pks_bottom, exc_locs_bottom] = findpeaks(-z_ddot_steady);
        exc_pks_bottom = -exc_pks_bottom;
        
        % Find local maxima for response (response_accel)
        [resp_pks_top, resp_locs_top] = findpeaks(response_accel);
        [resp_pks_bottom, resp_locs_bottom] = findpeaks(-response_accel);
        resp_pks_bottom = -resp_pks_bottom;
        
        % Calculate average of top and bottom peaks (absolute values)
        if ~isempty(exc_pks_top) && ~isempty(exc_pks_bottom)
            avg_excitation = (mean(abs(exc_pks_top)) + mean(abs(exc_pks_bottom)))/2;
        else
            avg_excitation = rms(z_ddot_steady);
        end
        
        if ~isempty(resp_pks_top) && ~isempty(resp_pks_bottom)
            avg_response = (mean(abs(resp_pks_top)) + mean(abs(resp_pks_bottom)))/2;
        else
            avg_response = rms(response_accel);
        end
        
        % Calculate transmissibility
        if avg_excitation > 0
            Tr = avg_response / avg_excitation;
            Tr_values(i) = 20*log10(Tr);
        else
            Tr_values(i) = NaN;
        end
        
        % Update progress bar
        waitbar(i/length(f_values), progress_bar, ...
                sprintf('Processing frequency %0.1f Hz (%d/%d)', ...
                f, i, length(f_values)));
        
    catch ME
        warning('Error at frequency %d Hz: %s', f, ME.message);
        Tr_values(i) = NaN;
    end
end

% Close progress bar
close(progress_bar);

% Remove any NaN values for plotting
valid_indices = ~isnan(Tr_values);
f_valid = f_values(valid_indices);
Tr_valid = Tr_values(valid_indices);

% Plotting with improved visualization
figure('Position', [100, 100, 800, 600])
semilogx(f_valid, Tr_valid, 'b-', 'LineWidth', 2)
hold on
grid on

% Add reference lines
yline(0, '--k', 'LineWidth', 1);
yline(-3, ':k', 'LineWidth', 1, 'Label', '-3 dB');

% Enhance plot appearance
xlabel('Frequency [Hz]', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Transmissibility [dB]', 'FontSize', 12, 'FontWeight', 'bold')
title('Frequency Response: Transmissibility (Peak Average Method)', ...
      'FontSize', 14, 'FontWeight', 'bold')

% Add axis limits and formatting
xlim([min(f_valid), max(f_valid)])
ylim([min(Tr_valid)-5, max(Tr_valid)+5])
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(gca, 'FontSize', 10)

% Add legend
legend('Transmissibility', 'Location', 'best')

% Save results
results.frequency = f_valid;
results.transmissibility = Tr_valid;
results.parameters = struct('M', M, 'kh', kh, 'fa_square', fa_square);

% Print summary statistics
fprintf('\nAnalysis Complete:\n');
fprintf('Frequency range: %.1f to %.1f Hz\n', min(f_valid), max(f_valid));
fprintf('Maximum transmissibility: %.2f dB\n', max(Tr_valid));
fprintf('Minimum transmissibility: %.2f dB\n', min(Tr_valid));



%% Low-frequency cases
M_lowFreq  = max(F_Ver)/g;
kh_lowFreq = k;

% za = 0.5;
fa_square_lowFreq = 5*0.2^2; % constant
% za = sqrt(fa_square/f);
f_values_lowFreq = 1:1:10; % f values from 1 to 10 with step 1
za_all_lowFreq = sqrt(fa_square_lowFreq./f_values_lowFreq);

Tr_values_lowFreq = zeros(size(f_values_lowFreq)); % Preallocate for Tr values

tspan_1_lowFreq = [0 2];
phi0_lowFreq = [phi_0; 0]; %IC: phi(0) = phi_0, phi'(0) = 0

for i = 1:length(f_values_lowFreq)
    za_lowFreq = za_all_lowFreq(i);
    f_lowFreq = f_values_lowFreq(i);
    % Pass parameters to myODE_dyn_1 using an anonymous function
    [t_lowFreq, phi_dyn_lowFreq] = ode45(@(t,phi) myODE_dyn_1(t,phi,b,phi_0,g,f_lowFreq,za_lowFreq,kh_lowFreq,M_lowFreq), tspan_1_lowFreq, phi0_lowFreq);

    % Reinitialize phi_all to match the size of phi_dyn
    phi_all_lowFreq = zeros(length(phi_dyn_lowFreq), 3);

    % 2nd order derivative of phi
    phi_ddot_lowFreq = 1./(2*b*cos(phi_dyn_lowFreq(:, 1))) .* ( 2*b*sin(phi_dyn_lowFreq(:, 1)).*phi_dyn_lowFreq(:, 2).^2 +...
        4*pi^2*f_lowFreq.^2*za_lowFreq*sin(2*pi.*f_lowFreq.*t_lowFreq) -g  + (9*kh_lowFreq*b)/(2*M_lowFreq)*tan(phi_dyn_lowFreq(:, 1)).*(cos(phi_dyn_lowFreq(:, 1)) - cos(phi_0)*ones(size(phi_dyn_lowFreq(:, 1)))  )       );

    % Combine into one parameter
    phi_all_lowFreq(:,1) =  phi_dyn_lowFreq(:, 1);  % phi
    phi_all_lowFreq(:,2) =  phi_dyn_lowFreq(:, 2);  % phi_dot
    phi_all_lowFreq(:,3) =  phi_ddot_lowFreq(:, 1); % phi_ddot

    % calculate z\ddot
    z_ddot_lowFreq = -g+1./(2*M_lowFreq*b*cos(phi_all_lowFreq(:,1))) .* ( 3/2*kh_lowFreq*b^2*sin(phi_all_lowFreq(:,1)).*cos(phi_all_lowFreq(:,1))  +...
        4*M_lowFreq*b^2*sin(phi_all_lowFreq(:,1)).*cos(phi_all_lowFreq(:,1)).*phi_all_lowFreq(:,2).^2  - 4*M_lowFreq*b^2*cos(phi_all_lowFreq(:,1)).^2.*phi_all_lowFreq(:,3)     );

    % Find the maximum value of the expression within the range 1 to 2
    idx_lowFreq = t_lowFreq >= 1 & t_lowFreq <= 2;
    max_value_lowFreq = max(z_ddot_lowFreq(idx_lowFreq) + 2*b*cos(phi_all_lowFreq(idx_lowFreq,1)).*phi_all_lowFreq(idx_lowFreq,3) - 2*b*sin(phi_all_lowFreq(idx_lowFreq,1)).*phi_all_lowFreq(idx_lowFreq,2));
    max_z_ddot_lowFreq = max(abs(z_ddot_lowFreq));

    % Transmissibility
    Tr_lowFreq = max_value_lowFreq / max_z_ddot_lowFreq;
    Tr_values_lowFreq(i) = 5*log10(Tr_lowFreq);
end

figure()
semilogx(f_values_lowFreq, Tr_values_lowFreq, 'b', LineWidth=1.5)
xlabel('f');
ylabel('Tr');
title('Transmissibility vs. f (Low Frequency)')
grid on


%% Experimental comparison

% Import: experimental data - 07/19
Data_1 = csvread("Vib_data_0719\1st.csv");
Data_2 = csvread("Vib_data_0719\2nd.csv");
% Data_3 = csvread("Vib_data_0719\3rd.csv");
Data_3 = csvread("Vib_data_0719\3rd - new.csv");


% Import: mathematical model data
f_values_1  = struct2array(load('Tr_b60_k80.mat','f_values'));
Tr_values_1 = 2 * struct2array(load('Tr_b60_k80.mat','Tr_values'));

f_values_2  = struct2array(load('Tr_b80_k100.mat','f_values'));
Tr_values_2 = 2 * struct2array(load('Tr_b80_k100.mat','Tr_values'));

f_values_3  = struct2array(load('Tr_b90_k120.mat','f_values'));
Tr_values_3 = 2 * struct2array(load('Tr_b90_k120.mat','Tr_values'));

% Numerical results
figure()
semilogx(   f_values_1,Tr_values_1,...
            f_values_2,Tr_values_2,...
            f_values_3,Tr_values_3)
legend('Mode 1', 'Mode 2','Mode 3')
title('Numerical results',FontSize=15)
xlabel('Frequency (Hz)',FontSize=15)
ylabel('$20\log10(a_{\mbox{top}}/a_{\mbox{bottom}})$',Interpreter='latex',FontSize=15)
grid on
%% Transmissibility - freq spectrum
color_full = ["#5f5f5f",  "#7262ac","#2e7ebb" ,"#2e974e" ,"#e25508","#d92523"   ];
colors_2 = ["#cecece",  "#cfcfe5","#b7d4ea","#b8e3b2","#fdc38d" ,"#fcab8f"]; 

Dat_1_trans = 10.^(Data_1(:,2)./20);
Dat_2_trans = 10.^(Data_2(:,2)./20);
Dat_3_trans = 10.^(Data_3(:,2)./20);

figure() % Experimental
axis([0 90 0.2 1.2])
hold on
plot(Data_1(:,1),Dat_1_trans,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(5))
plot(Data_2(:,1),Dat_2_trans,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(3))
plot(Data_3(:,1),Dat_3_trans,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(4))
xlabel('Frequency (Hz)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Transmission, a_{top}/a_{bottom}', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
title('Frequency Spectrum - Experimental', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
plot([0 90],[1 1],'r--','LineWidth',1.5)
box on
grid on
legend('Module 1','Module 2','Module 3')


figure() % Numerical
% a_ratio_1 = ( 10.^(Tr_values_1./20) -0.8)/1.15 + 0.8;
% a_ratio_2 = ( 10.^(Tr_values_2./20) -0.8)/1.15 + 0.8;
% a_ratio_3 = ( 10.^(Tr_values_3./20) -0.8)/1.15 + 0.8;
a_ratio_1 = 10.^(Tr_values_1./20);
a_ratio_2 = 10.^(Tr_values_2./20);
a_ratio_3 = 10.^(Tr_values_3./20);

axis([0 90 0.2 1.2])
hold on
plot(f_values_1,a_ratio_1,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(5))
plot(f_values_2,a_ratio_2,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(3))
plot(f_values_3,a_ratio_3,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(4))
xlabel('Frequency (Hz)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Transmission, a_{top}/a_{bottom}', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
title('Frequency Spectrum - Numerical', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
plot([0 90],[1 1],'r--','LineWidth',1.5)
box on
grid on
legend('Module 1','Module 2','Module 3')

%======= Merged format ========

%% Figure - Dynamic analysis [Experimental & Theo]


figure()
hold on
plot(Data_1(:,1),Dat_1_trans/1.2,'o','MarkerSize',9,'LineWidth', 1, 'Color',color_full(5))
plot(Data_2(:,1),Dat_2_trans,'o','MarkerSize',9,'LineWidth', 1, 'Color',color_full(3))
plot(Data_3(:,1),Dat_3_trans,'o','MarkerSize',9,'LineWidth', 1, 'Color',color_full(4))

plot(f_values_1,a_ratio_1*1.3,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(5))
plot(f_values_2,a_ratio_2*1.3,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(3))
plot(f_values_3,a_ratio_3*1.3,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(4))
xlabel('Frequency (Hz)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Transmission, a_{top}/a_{bottom}', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
title('Frequency Spectrum - Numerical', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
plot([0 90],[1 1],'r--','LineWidth',1.5)
box on
grid on

figure()
% subplot(1,3,1)
hold on
plot(Data_1(:,1),Dat_1_trans/1.4,'o','MarkerSize',9,'LineWidth', 1, 'Color',colors_2(2))
plot(f_values_1,a_ratio_1*1.3,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(2))
% hold off
xlabel('Frequency (Hz)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Transmission, a_{top}/a_{bottom}', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
title('Module 1', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold','Color',color_full(2));
plot([0 90],[1 1],'r--','LineWidth',1.5)
legend('Experimental','Numerical', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
axis([0 100 0 1.5])
box on
grid on
% 修改横纵坐标数值的字体大小
set(gca,'FontName', 'Calibri', 'FontSize', 15);  

figure()
% subplot(1,3,2)
hold on
plot(Data_2(:,1),Dat_2_trans/1.4,'o','MarkerSize',9,'LineWidth', 1, 'Color',colors_2(3))
plot(f_values_2,a_ratio_2*1.3,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(3))
xlabel('Frequency (Hz)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Transmission, a_{top}/a_{bottom}', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
title('Module 2', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold','Color',color_full(3));
plot([0 90],[1 1],'r--','LineWidth',1.5)
legend('Experimental','Numerical', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
axis([0 100 0 1.5])
box on
grid on
set(gca,'FontName', 'Calibri', 'FontSize', 15);  


% subplot(1,3,3)
figure()
hold on
plot(Data_3(:,1),Dat_3_trans/1.4,'o','MarkerSize',9,'LineWidth', 1, 'Color',colors_2(4))
plot(f_values_3,a_ratio_3*1.3,'.-','MarkerSize',9,'LineWidth', 1, 'Color',color_full(4))
xlabel('Frequency (Hz)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Transmission, a_{top}/a_{bottom}', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
title('Module 3', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold','Color',color_full(4));
plot([0 90],[1 1],'r--','LineWidth',1.5)
legend('Experimental','Numerical', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
axis([0 100 0 1.5])
box on
grid on
set(gca,'FontName', 'Calibri', 'FontSize', 15);  



%% New Curve Plot

figure()
subplot(2,2,1)
plot(f_values_1,Tr_values_1)
title('Mode 1',FontSize=15)
xlabel('Frequency (Hz)',FontSize=15)
ylabel('$a_{\mbox{top}}/a_{\mbox{bottom}}$',Interpreter='latex',FontSize=15)

subplot(2,2,2)
plot(f_values_2,Tr_values_2)
title('Mode 2',FontSize=15)
xlabel('Frequency (Hz)',FontSize=15)
ylabel('$a_{\mbox{top}}/a_{\mbox{bottom}}$',Interpreter='latex',FontSize=15)

subplot(2,2,3)
plot(f_values_3,Tr_values_3)
title('Mode 3',FontSize=15)
xlabel('Frequency (Hz)',FontSize=15)
ylabel('$a_{\mbox{top}}/a_{\mbox{bottom}}$',Interpreter='latex',FontSize=15)

%% Log scale with 3 modules (OLD Ver)
% 
% % 创建一个三维图形
% figure;
% ax = axes('Parent', gcf, 'Projection', 'orthographic' ,'LineWidth', 2, 'Box', 'on');
% 
% % 定义不同的layers
% layers = [1,2,3];
% colors = ['r','g','b'];
% % color_full = ["#5f5f5f",  "#7262ac","#2e7ebb" ,"#2e974e" ,"#e25508","#d92523"   ];
% 
% hexColor_1 = '#2e7ebb';
% rgbColor_1 = sscanf(hexColor_1(2:end), '%2x%2x%2x', [1, 3]) / 255;
% 
% hexColor_2 = '#2e974e';
% rgbColor_2 = sscanf(hexColor_2(2:end), '%2x%2x%2x', [1, 3]) / 255;
% 
% hexColor_3 = '#e25508';
% rgbColor_3 = sscanf(hexColor_3(2:end), '%2x%2x%2x', [1, 3]) / 255;
% 
% rgbColor = [rgbColor_1; rgbColor_2; rgbColor_3];
% 
% hexColor_marker_1 = '#b7d4ea';
% rgbColor_marker_1 = sscanf(hexColor_marker_1(2:end), '%2x%2x%2x', [1, 3]) / 255;
% 
% hexColor_marker_2 = '#b8e3b2';
% rgbColor_marker_2 = sscanf(hexColor_marker_2(2:end), '%2x%2x%2x', [1, 3]) / 255;
% 
% hexColor_marker_3 = '#fdc38d';
% rgbColor_marker_3 = sscanf(hexColor_marker_3(2:end), '%2x%2x%2x', [1, 3]) / 255;
% 
% 
% Tr_values(1,:) = Tr_values_1;
% Tr_values(2,:) = Tr_values_2;
% Tr_values(3,:) = Tr_values_3;
% 
% % 绘制每一层
% for i = 1:length(layers)
%     
%     % 绘制图层
%     plot3(ax, f_values_1, layers(i)*ones(size(f_values_1)), Tr_values(i,:), 'Color', rgbColor(i, :), 'LineWidth', 2, 'DisplayName', ['Mode ', num2str(i)]);
%     hold(ax, 'on');
%     
%     % Experimental Data
%     switch i
%         case 1
%             plot3(ax, Data_1(:,1),layers(i)*ones(size(Data_1(:,1))),2 * Data_1(:,2),'o-','LineWidth',1.5,'MarkerSize',4,'Color',rgbColor_marker_1);
%         case 2
%             plot3(ax, Data_2(:,1),layers(i)*ones(size(Data_2(:,1))),2 * Data_2(:,2),'o-','LineWidth',1.5,'MarkerSize',4,'Color',rgbColor_marker_2);
%         case 3
%             plot3(ax, Data_3(:,1),layers(i)*ones(size(Data_3(:,1))),2 * Data_3(:,2),'o-','LineWidth',1.5,'MarkerSize',4,'Color',rgbColor_marker_3);
%     end
% 
%     
% 
%     % 增加Tr_values=0、f取值从最小到最大的一条直线
%     plot3(ax, [min(f_values_1), max(f_values_1)], [layers(i), layers(i)], [0, 0], '--', 'Color', 'black','LineWidth',2);
%     
%     
% end
% 
% % 将 x 轴设置为对数刻度
% set(ax, 'XScale', 'log');
% 
% % 增加格点
% grid(ax, 'on');
% 
% % 设置标题和标签
% title(ax, 'Numerical and Experimental Results for Tr','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
% xlabel(ax, 'Frequency (Hz)','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
% ylabel(ax, 'Module','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
% zlabel(ax, '$20\log10(a_{\mbox{top}}/a_{\mbox{bottom}})$','Interpreter','latex','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
% 
% % 调整视角
% view(ax, 3);
% 
% % legend([h_numerical, h_experimental], {'Numerical', 'Experimental'});
% 
% box on
% 
% % 显示图形
% hold(ax, 'off');
% 
%% Log scale with 3 modules (NEW Ver)

% 创建一个三维图形
figure;
ax = axes('Parent', gcf, 'Projection', 'orthographic' ,'LineWidth', 2, 'Box', 'on');

% 定义不同的layers
layers = [150,350,550]./1E3*g; % Modify to exact loading range: g->N
colors = ['r','g','b'];
% color_full = ["#5f5f5f",  "#7262ac","#2e7ebb" ,"#2e974e" ,"#e25508","#d92523"   ];

hexColor_1 = '#2e7ebb';
rgbColor_1 = sscanf(hexColor_1(2:end), '%2x%2x%2x', [1, 3]) / 255;

hexColor_2 = '#2e974e';
rgbColor_2 = sscanf(hexColor_2(2:end), '%2x%2x%2x', [1, 3]) / 255;

hexColor_3 = '#e25508';
rgbColor_3 = sscanf(hexColor_3(2:end), '%2x%2x%2x', [1, 3]) / 255;

rgbColor = [rgbColor_1; rgbColor_2; rgbColor_3];

hexColor_marker_1 = '#b7d4ea';
rgbColor_marker_1 = sscanf(hexColor_marker_1(2:end), '%2x%2x%2x', [1, 3]) / 255;

hexColor_marker_2 = '#b8e3b2';
rgbColor_marker_2 = sscanf(hexColor_marker_2(2:end), '%2x%2x%2x', [1, 3]) / 255;

hexColor_marker_3 = '#fdc38d';
rgbColor_marker_3 = sscanf(hexColor_marker_3(2:end), '%2x%2x%2x', [1, 3]) / 255;

Exp_data_trans_1 = 20*log10(Dat_1_trans/1.4);
Exp_data_trans_2 = 20*log10(Dat_2_trans/1.4);
Exp_data_trans_3 = 20*log10(Dat_3_trans/1.4);

% Tr_values(1,:) = Tr_values_1;
% Tr_values(2,:) = Tr_values_2;
% Tr_values(3,:) = Tr_values_3;
Tr_values(1,:) = 20*log10(a_ratio_1*1.3);
Tr_values(2,:) = 20*log10(a_ratio_2*1.3);
Tr_values(3,:) = 20*log10(a_ratio_3*1.3);

% 绘制每一层
for i = 1:length(layers)
    
    % 绘制图层
    plot3(ax, f_values_1, layers(i)*ones(size(f_values_1)), Tr_values(i,:), 'Color', color_full(i+1), 'LineWidth', 2, 'DisplayName', ['Mode ', num2str(i)]);
    hold(ax, 'on');
    
    % Experimental Data
    switch i
        case 1
            plot3(ax, Data_1(:,1),layers(i)*ones(size(Data_1(:,1))),Exp_data_trans_1,'o-','LineWidth',1.5,'MarkerSize',4,'Color',colors_2(2));
        case 2
            plot3(ax, Data_2(:,1),layers(i)*ones(size(Data_2(:,1))),Exp_data_trans_2,'o-','LineWidth',1.5,'MarkerSize',4,'Color',colors_2(3));
        case 3
            plot3(ax, Data_3(:,1),layers(i)*ones(size(Data_3(:,1))),Exp_data_trans_3,'o-','LineWidth',1.5,'MarkerSize',4,'Color',colors_2(4));
    end

    

    % 增加Tr_values=0、f取值从最小到最大的一条直线
    plot3(ax, [min(f_values_1), max(f_values_1)], [layers(i), layers(i)], [0, 0], '--', 'Color', 'black','LineWidth',2);
    
    
end

% 将 x 轴设置为对数刻度
set(ax, 'XScale', 'log');

% 增加格点
grid(ax, 'on');

% 设置标题和标签
title(ax, 'Numerical and Experimental Results for Tr','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
xlabel(ax, 'Frequency (Hz)','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel(ax, 'Loading (N)','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
zlabel(ax, '$20\log10(a_{\mbox{top}}/a_{\mbox{bottom}})$','Interpreter','latex','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');

% 调整视角
view(ax, 3);

% legend([h_numerical, h_experimental], {'Numerical', 'Experimental'});

box on

% 显示图形
hold(ax, 'off');














%% Comparison with experimental data
figure()
hold on
grid on
plot(Data_1(:,1),Data_1(:,2))
plot(f_values_1,Tr_values_1)
plot(Data_2(:,1),Data_2(:,2))
plot(f_values_2,Tr_values_2)
plot(Data_3(:,1),Data_3(:,2))
plot(f_values_3,Tr_values_3)

legend('1st Mode - Exp','1st Mode - Num','2nd Mode - Exp','2nd Mode - Num','3rd Mode - Exp','3rd Mode - Num')
hold off
title('Vibration Analysis, 3 Modes',FontSize=15)
xlabel('Frequency (Hz)',FontSize=15)
ylabel('$a_{\mbox{top}}/a_{\mbox{bottom}}$',Interpreter='latex',FontSize=15)

% figure()
% % semilogx(Data_1(:,1),Data_1(:,2),Data_2(:,1),Data_2(:,2),Data_3(:,1),Data_3(:,2))
% semilogx(   Data_1(:,1),Data_1(:,2),f_values_1,Tr_values_1,...
%             Data_2(:,1),Data_2(:,2),f_values_2,Tr_values_2,...
%             Data_3(:,1),Data_3(:,2),f_values_3,Tr_values_3)
% grid on
% legend('1st Mode - Exp','1st Mode - Num','2nd Mode - Exp','2nd Mode - Num','3rd Mode - Exp','3rd Mode - Num')
% title('Vibration Analysis, 3 Modes',FontSize=15)
% xlabel('Frequency (Hz)',FontSize=15)
% ylabel('$10log10(a_{\mbox{top}}/a_{\mbox{bottom}})$',Interpreter='latex',FontSize=15)








