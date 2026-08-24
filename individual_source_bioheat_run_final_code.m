clc; clear all; close all;
% Read the following instructions carefully before running the script. 
% This script is designed to perform a multi-wavelength Monte Carlo simulation of light transport in tissue, followed by a bio-heat transfer analysis using the Pennes bioheat equation. The results will be saved in an Excel file and visualized in several figures.
% Whoever runs this script, please ensure that the MCXLAB engine is installed and properly configured in your MATLAB path. 
% The script will automatically check for the presence of MCXLAB and will throw an error if it is not found.
% Please make a copy of the standalone spectrum CSV file (e.g., TLS_1450_Ref-1.csv) in the same directory as this script or provide the full path to the file.
% Make sure you have the necessary permissions to write files in the directory where this script is located, as it will create a results folder and save output files there.
% copy this script to a new folder and run it from there to avoid any path conflicts with other scripts or functions.
% any changes to the tissue optical properties or bio-thermal constants should be made in the respective sections of the script, and not in the helper function, to ensure consistency across simulations.
% If you change the mxc outut type to fluence rate from fluence, please ensure that the units are consistent throughout the script, 
% especially while multiplying by spectral power you have to multiply by energy in Joules and not in Watts, as the fluence rate is in W/mm^2 and the fluence is in J/mm^2 and fluence rate output is in 1/mm^2 
% mcx output for fluence in 1/mm^2 and for fluence rate it is 1/mm^2 s.
% make sure the boundary conditions in the FDM solver are consistent with the physical scenario you are simulating, as incorrect boundary conditions can lead to non-physical results.
% If you read above and understand the above instructions, you can proceed to run the script. 
% The script will perform a multi-wavelength Monte Carlo simulation of light transport in tissue, followed by a bio-heat transfer analysis using the Pennes bioheat equation. The results will be saved in an Excel file and visualized in several figures.
% Now you are good to go. Enjoy the simulation and analysis!

% =============================================================================
% 0. AUTOMATIC PATH RESOLUTION & MCXLAB CHECK
% =============================================================================
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
parent_dir = fileparts(script_dir);

addpath(genpath(script_dir));
if exist(parent_dir, 'dir'), addpath(genpath(parent_dir)); end

mcx_explicit_path = 'C:\Users\raj.thakur\OneDrive - University College Cork\RI Pathway NASCENT - Documents\General\09_Data\Experiment\Bio-Heat Simulation\working_programs\MCXLAB\mcx-master\mcxlab';
if exist(mcx_explicit_path, 'dir'), addpath(genpath(mcx_explicit_path)); end

if exist('mcxlab', 'file') == 2 || exist('mcxlab', 'file') == 3
    fprintf('[SUCCESS] MCXLAB engine detected.\n');
else
    error('MCXLAB engine not found. Verify installation directory.');
end

% =============================================================================
% 1. GRID & SPATIAL DEFINITION
% =============================================================================
tissue_x = 10;   % mm
tissue_y = 10;   % mm
tissue_z = 20;   % mm
dx = 0.1;        % Grid step x (mm)
dy = 0.1;        % Grid step y (mm)
dz = 0.2;        % Grid step z (mm)

nx = round(tissue_x / dx);
ny = round(tissue_y / dy);
nz = round(tissue_z / dz);

xx = linspace(0.0, tissue_x*1e-3, nx); 
yy = linspace(0.0, tissue_y*1e-3, ny);
zz = linspace(0.0, tissue_z*1e-3, nz);

dx_m = dx * 1e-3; % meters
dy_m = dy * 1e-3; % meters
dz_m = dz * 1e-3; % meters

cx = round(nx / 2); % Beam center index x
cy = round(ny / 2); % Beam center index y

% =============================================================================
% 2. BIO-THERMAL CONSTANTS
% =============================================================================
k       = 0.4975;   % Thermal conductivity (W/(m*K))
c_b     = 3664;     % Blood specific heat (J/(kg*K))
rho_b   = 1100;     % Blood density (kg/m^3)
w_b     = 1.58e-3;  % Blood perfusion rate (1/s)
Q_m     = 678;      % Metabolic heat generation (W/m^3)
T_b     = 37.0;     % Arterial blood baseline temperature (°C)
T_amb   = 35.0;     % Ambient temperature (°C)
h       = 4.5;      % Heat transfer coefficient (W/(m^2*K))
epsilon = 0.9;      % Emissivity
sigma   = 5.67e-8;  % Stefan-Boltzmann constant (W/(m^2*K^4))

S_perfusion = rho_b * c_b * w_b; % Perfusion damping coefficient

% =============================================================================
% 3. STANDALONE SPECTRUM FILE LOADING & LOSS COMPENSATION
% =============================================================================
% standalone_filename = 'Lamp_Ref-1.csv';
% standalone_filename = 'LED-1450_Ref-1.csv';
% standalone_filename = 'TDFA19002_Ref-1.csv';
standalone_filename = 'TLS_1450_Ref-1.csv';
% standalone_filename = 'TLS_1470_Ref-1.csv';

bin_step_nm = 5;  % Spectral binning resolution (nm)
loss_dB     = 3.75; % Optical loss compensation (dB)

csv_fullpath = which(standalone_filename);
if isempty(csv_fullpath)
    files = dir(fullfile(parent_dir, '**', standalone_filename));
    if ~isempty(files)
        csv_fullpath = fullfile(files(1).folder, files(1).name);
    else
        csv_fullpath = standalone_filename;
    end
end

fprintf('--> Loading Standalone Spectrum File: %s\n', csv_fullpath);
data_table = readtable(csv_fullpath);

raw_wavelengths = double(data_table{:, 1})';
raw_power_dBm   = double(data_table{:, 2})';

if max(raw_wavelengths) < 1.0
    raw_wavelengths = raw_wavelengths * 1e9;
    fprintf('--> Wavelengths detected in METERS. Converted to nm.\n');
end

% Compensate for loss and convert to linear Watts
raw_power_dBm_compensated = raw_power_dBm + loss_dB ;
raw_power_W = 10.^(raw_power_dBm_compensated / 10) * 1e-3;
raw_power_W(raw_power_dBm <= -120 | isnan(raw_power_W)) = 0;
source_label = strrep(standalone_filename, '.csv', '');

% --- SPECTRAL FILTER & PREALLOCATED BINNING ---
active_mask = raw_power_W > 1e-12;
active_wl   = raw_wavelengths(active_mask);
active_pw   = raw_power_W(active_mask);

bin_edges = floor(min(active_wl)) : bin_step_nm : ceil(max(active_wl));
num_bins  = length(bin_edges) - 1;

binned_wl_temp = zeros(1, num_bins);
binned_pw_temp = zeros(1, num_bins);
count = 0;

for b = 1:num_bins
    band_mask = (active_wl >= bin_edges(b)) & (active_wl < bin_edges(b+1));
    if any(band_mask)
        count = count + 1;
        binned_wl_temp(count) = mean(active_wl(band_mask));
        binned_pw_temp(count) = sum(active_pw(band_mask));
    end
end

binned_wl = binned_wl_temp(1:count);
binned_pw = binned_pw_temp(1:count);
num_bands = length(binned_wl);

fprintf('\n>>> Target Source: %s\n', source_label);
fprintf('>>> Total Integrated Power (Compensated): %.6f mW (%.2f uW) across %d spectral bands\n', ...
        sum(binned_pw)*1e3, sum(binned_pw)*1e6, num_bands);

% =============================================================================
% 4. MULTI-WAVELENGTH MCX SIMULATION & HEAT ACCUMULATION
% =============================================================================
Q_laser = zeros(nx, ny, nz);
total_fluence = zeros(nx, ny, nz); % Total physical fluence rate (W/mm^2)
beam_radius_mm = 0.5;

delta_T_spectral    = zeros(1, num_bands);
T_spectral          = zeros(1, num_bands);
T_diffused_spectral = zeros(1, num_bands);
mu_a_record         = zeros(1, num_bands);
mu_s_record         = zeros(1, num_bands);

fprintf('\n=== Launching GPU Monte Carlo Light Transport (MCXLAB) ===\n');
for i = 1:num_bands
    wl_nm = binned_wl(i);
    P_W   = binned_pw(i);
    
    [mu_a, mu_s, g, n] = get_tissue_optical_properties(wl_nm);
    mu_a_record(i) = mu_a;
    mu_s_record(i) = mu_s;
    
    if mu_a == 0 || P_W == 0
        delta_T_spectral(i)    = 0.0;
        T_spectral(i)          = T_b;
        T_diffused_spectral(i) = T_b;
        continue;
    end
    
    cfg.nphoton  = 1e6;
    cfg.vol      = uint8(ones(nx, ny, nz));
    cfg.srctype  = 'gaussian';
    cfg.srcpos   = [nx/2, ny/2, 1];
    cfg.srcdir   = [0, 0, 1];
    cfg.unitinmm = dx;
    cfg.srcparam1= [beam_radius_mm/dx, 0, 0, 0];
    cfg.prop     = [0, 0, 1, 1; mu_a, mu_s, g, n];
    cfg.tstart   = 0;
    cfg.tend     = 5e-9;
    cfg.tstep    = cfg.tend;
    cfg.gpuid    = 1;
    cfg.outputtype = 'fluence';  
    
    fluence = mcxlab(cfg);
    fluence_mcx = double(fluence.data);

    % Accumulate total physical optical fluence (W/mm^2)
    total_fluence = total_fluence + (fluence_mcx .* P_W);
    
    Q_lambda = (mu_a .* fluence_mcx .* P_W) * 1e9; % Volumetric heat (W/m^3)
    Q_laser  = Q_laser + Q_lambda;
    
    % Un-diffused local peak
    delta_T_spectral(i) = max(Q_lambda(:)) / S_perfusion;
    T_spectral(i)       = T_b + delta_T_spectral(i);
    
    % Single-band 3D PDE-diffused peak
    T_band = T_b * ones(nx, ny, nz);
    denom_b = 2/dx_m^2 + 2/dy_m^2 + 2/dz_m^2 + S_perfusion/k;
    err_b = Inf;
    
    while err_b > 1e-3
        T_old_b = T_band;
        lap_b = (T_band(3:end, 2:end-1, 2:end-1) + T_band(1:end-2, 2:end-1, 2:end-1)) / dx_m^2 + ...
                (T_band(2:end-1, 3:end, 2:end-1) + T_band(2:end-1, 1:end-2, 2:end-1)) / dy_m^2 + ...
                (T_band(2:end-1, 2:end-1, 3:end) + T_band(2:end-1, 2:end-1, 1:end-2)) / dz_m^2;
        
        T_band(2:end-1, 2:end-1, 2:end-1) = ...
            (lap_b + (S_perfusion/k)*T_b + (Q_m + Q_lambda(2:end-1, 2:end-1, 2:end-1))/k) / denom_b;
        
        T_band(:, :, end) = T_b;
        T_band(1, :, :)   = T_band(2, :, :);   T_band(end, :, :) = T_band(end-1, :, :);
        T_band(:, 1, :)   = T_band(:, 2, :);   T_band(:, end, :) = T_band(:, end-1, :);
        T_band(:, :, 1)   = T_band(:, :, 2);
        
        err_b = max(abs(T_band(:) - T_old_b(:)));
    end
    T_diffused_spectral(i) = max(T_band(:));
    
    fprintf('[Band %2d/%d] λ = %.1f nm | P = %.6f uW | mu_a = %.2f cm^-1 | ΔT = %.6f °C\n', ...
            i, num_bands, wl_nm, P_W*1e6, mu_a*10, delta_T_spectral(i));
end

%====================================================================================

% =============================================================================
% 5. 3D PENNES BIOHEAT SOLVER (PURE EQUILIBRIUM CONVERGENCE)
% =============================================================================
T = T_b * ones(nx, ny, nz);
denom = 2/dx_m^2 + 2/dy_m^2 + 2/dz_m^2 + S_perfusion/k;
h_eff = h + 4 * epsilon * sigma * (T_amb)^3; % Linearized Robin+Radiation
fprintf('\n--- Solving 3D Pennes Bioheat Equation to Equilibrium ---\n');
err = Inf; 
it = 0; 
tol = 1e-6; % Equilibrium criterion: delta_T < 10^-6 °C (1 micro-Kelvin)

% Dynamic tracking arrays
iter_history = [];
time_history = [];
temp_history = [];
err_history  = [];

t_solver_start = tic;

% LOOP RUNS PURELY UNTIL THERMAL EQUILIBRIUM (No artificial iteration cap)
while err > tol
    it = it + 1;
    T_old = T;
    
    % 3D Laplacian (Conduction term)
    lap = (T(3:end, 2:end-1, 2:end-1) + T(1:end-2, 2:end-1, 2:end-1)) / dx_m^2 + ...
          (T(2:end-1, 3:end, 2:end-1) + T(2:end-1, 1:end-2, 2:end-1)) / dy_m^2 + ...
          (T(2:end-1, 2:end-1, 3:end) + T(2:end-1, 2:end-1, 1:end-2)) / dz_m^2;
    
    % Gauss-Seidel / Jacobi update
    T(2:end-1, 2:end-1, 2:end-1) = ...
        (lap + (S_perfusion/k)*T_b + (Q_m + Q_laser(2:end-1, 2:end-1, 2:end-1))/k) / denom;
    
    % Boundary Conditions (Dirichlet at bottom, Neumann adiabatic on sides/top)
    T(:, :, end) = T_b;
    T(1, :, :)   = T(2, :, :);   T(end, :, :) = T(end-1, :, :);
    T(:, 1, :)   = T(:, 2, :);   T(:, end, :) = T(:, end-1, :);
    T(:, :, 1)   = T(:, :, 2);

    % Top Surface (z = 1) Convective + Radiative Cooling
    T(:, :, 1)   = (k * T(:, :, 2) + dz_m * h_eff * T_amb) / (k + dz_m * h_eff);
    
    % Maximum pointwise temperature variation across the 3D grid
    err = max(abs(T(:) - T_old(:)));
    
    % Record convergence trajectory
    iter_history(it) = it;
    time_history(it) = toc(t_solver_start);
    temp_history(it) = max(T(:));
    err_history(it)  = err;
    
    if mod(it, 250) == 0
        fprintf('Iteration %4d (t = %5.2f s): Peak Temp = %.6f °C | Residual = %.3e °C\n', ...
                it, time_history(it), temp_history(it), err);
    end
end
elapsed_sec = toc(t_solver_start);

% =============================================================================
% VERIFY GLOBAL ENERGY CONSERVATION AT EQUILIBRIUM
% =============================================================================
voxel_vol_m3 = dx_m * dy_m * dz_m;
total_heat_generated = sum(Q_laser(:) + Q_m) * voxel_vol_m3; % Watts
total_heat_removed   = sum(S_perfusion * (T(:) - T_b)) * voxel_vol_m3; % Watts
net_energy_imbalance = abs(total_heat_generated - total_heat_removed);

final_peak_T  = max(T(:));
final_delta_T = final_peak_T - T_b;

fprintf('\n=================================================================\n');
fprintf(' [EQUILIBRIUM REACHED AUTONOMOUSLY]\n');
fprintf('  * Total Iterations to Equilibrium : %d\n', it);
fprintf('  * Total Solver Runtime            : %.2f seconds\n', elapsed_sec);
fprintf('  * Peak Steady-State Temperature   : %.6f °C (ΔT = +%.6f °C)\n', final_peak_T, final_delta_T);
fprintf('  * Pointwise Convergence Residual  : %.2e °C (Threshold: %.0e °C)\n', err, tol);
fprintf('  * Net Energy Balance Residual     : %.4e Watts\n', net_energy_imbalance);
fprintf('=================================================================\n');

% =============================================================================
% 6. DIRECTORY CREATION & STRUCTURED EXCEL EXPORT (NON-OVERWRITING)
% =============================================================================
timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
results_base_dir = fullfile(script_dir, 'Simulation_Results');
if ~exist(results_base_dir, 'dir'), mkdir(results_base_dir); end

run_folder = fullfile(results_base_dir, sprintf('%s_%s', source_label, timestamp_str));
if ~exist(run_folder, 'dir'), mkdir(run_folder); end

excel_filename = fullfile(run_folder, sprintf('%s_BioHeat_Results_%s.xlsx', source_label, timestamp_str));
fprintf('\n--> Saving all structured simulation data to Excel:\n    %s\n', excel_filename);

% --- Sheet 1: Metadata Summary Table ---
summary_names = { ...
    'Source_Label'; ...
    'Input_File'; ...
    'Date_Executed'; ...
    'Total_Power_mW'; ...
    'Total_Power_uW'; ...
    'Loss_Compensated_dB'; ...
    'Active_Spectral_Bands'; ...
    'Tissue_Dimensions_mm'; ...
    'Voxel_Size_mm'; ...
    'Peak_Steady_State_Temp_C'; ...
    'Max_Temp_Elevation_DeltaT_C'; ...
    'Convergence_Iterations'; ...
    'Solver_Runtime_Seconds' ...
};

summary_values = { ...
    source_label; ...
    standalone_filename; ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'); ...
    sum(binned_pw)*1e3; ...
    sum(binned_pw)*1e6; ...
    loss_dB; ...
    num_bands; ...
    sprintf('%d x %d x %d', tissue_x, tissue_y, tissue_z); ...
    sprintf('%.2f x %.2f x %.2f', dx, dy, dz); ...
    final_peak_T; ...
    final_delta_T; ...
    it; ...
    elapsed_sec ...
};
T_summary = table(summary_names, summary_values, 'VariableNames', {'Parameter', 'Value'});
writetable(T_summary, excel_filename, 'Sheet', 'Summary');

% --- Sheet 2: Spectral Results Table ---
T_spectral_table = table( ...
    binned_wl', ...
    binned_pw', ...
    (binned_pw*1e3)', ...
    (binned_pw*1e6)', ...
    mu_a_record', ...
    mu_s_record', ...
    T_spectral', ...
    T_diffused_spectral', ...
    'VariableNames', { ...
        'Wavelength_nm', ...
        'Power_Watts', ...
        'Power_mW', ...
        'Power_uW', ...
        'mu_a_mm_inv', ...
        'mu_s_mm_inv', ...
        'UnDiffused_Local_Temp_C', ...
        'Diffused_3D_PDE_Temp_C' ...
    });
writetable(T_spectral_table, excel_filename, 'Sheet', 'Spectral_Data');

% --- Sheet 3: Axial Depth Profile Table (Including Fluence) ---
z_mm = zz * 1e3;
T_axial = squeeze(T(cx, cy, :));
Q_axial = squeeze(Q_laser(cx, cy, :));
Phi_axial = squeeze(total_fluence(cx, cy, :));

T_depth_table = table( ...
    z_mm', ...
    T_axial, ...
    Q_axial, ...
    Q_axial * 1e-3, ...
    Phi_axial, ...
    log10(max(Phi_axial, 1e-15)), ...
    'VariableNames', { ...
        'Depth_z_mm', ...
        'Temperature_C', ...
        'Absorbed_Heat_Q_W_per_m3', ...
        'Absorbed_Heat_Q_kW_per_m3', ...
        'Fluence_Rate_W_per_mm2', ...
        'log10_Fluence_Rate' ...
    });
writetable(T_depth_table, excel_filename, 'Sheet', 'Depth_Profile');


% --- Sheet 4: FDM Convergence & Time History Table ---
T_conv_table = table( ...
    iter_history', ...
    time_history', ...
    temp_history', ...
    err_history', ...
    'VariableNames', { ...
        'Iteration', ...
        'Elapsed_Time_Seconds', ...
        'Peak_Temperature_C', ...
        'Pointwise_Error_Residual_C' ...
    });
writetable(T_conv_table, excel_filename, 'Sheet', 'FDM_Convergence');

fprintf('[SUCCESS] Multi-tab Excel workbook written successfully!\n');

% =============================================================================
% 7. VISUALIZATION DASHBOARDS & AUTOMATIC IMAGE SAVING
% =============================================================================
cx = round(nx / 2); cy = round(ny / 2);
x_mm = xx * 1e3; y_mm = yy * 1e3;

% --- Figure 1: Combined Bio-Thermal Analysis Dashboard ---
fig1 = figure('Name', sprintf('Combined Bio-Thermal Analysis [%s]', source_label), ...
              'Color', 'w', 'Position', [50, 50, 1100, 750]);

subplot(2, 2, 1);
plot(binned_wl, binned_pw * 1e6, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on; xlabel('Wavelength \lambda (nm)', 'FontSize', 11); ylabel('Optical Power (\muW)', 'FontSize', 11);
title(sprintf('Compensated Power Spectrum [%s]', source_label), 'FontSize', 12);

subplot(2, 2, 2);
plot(binned_wl, T_spectral, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on; xlabel('Wavelength \lambda (nm)', 'FontSize', 11); ylabel('Absolute Temperature T (°C)', 'FontSize', 11);
title('Absolute Temperature T vs Wavelength \lambda', 'FontSize', 12);

subplot(2, 2, 3);
plot(z_mm, T_axial, 'r-', 'LineWidth', 2.5);
grid on; xlabel('Tissue Depth z (mm)', 'FontSize', 11); ylabel('Temperature T (°C)', 'FontSize', 11);
title(sprintf('Axial Temp Profile (Peak T = %.6f °C)', final_peak_T), 'FontSize', 12);

subplot(2, 2, 4);
imagesc(x_mm, z_mm, squeeze(T(:, cy, :))');
axis image; cb = colorbar; ylabel(cb, 'Temperature (°C)', 'FontSize', 11);
colormap(jet); xlabel('Lateral Distance x (mm)', 'FontSize', 11); ylabel('Depth z (mm)', 'FontSize', 11);
title(sprintf('2D XZ Thermal Cross-Section - %s', source_label), 'FontSize', 12);

sgtitle(sprintf('Integrated Multi-Wavelength Bio-Thermal Analysis: %s', source_label), ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(run_folder, sprintf('%s_Dashboard_%s.png', source_label, timestamp_str)));

% --- Figure 2: Lateral Diffusion Profiles ---
[~, cz] = max(T_axial);
z_peak_mm = z_mm(cz);

fig2 = figure('Name', 'Lateral Thermal Diffusion Profiles', 'Color', 'w', 'Position', [80, 80, 1000, 450]);
subplot(1, 2, 1);
plot(x_mm, squeeze(T(:, cy, cz)), 'r-', 'LineWidth', 2.5);
grid on; xlabel('Lateral Distance x (mm)', 'FontSize', 11); ylabel('Temperature T (°C)', 'FontSize', 11);
title(sprintf('Lateral Profile T(x) at y = %.1f mm, z = %.2f mm', y_mm(cy), z_peak_mm), 'FontSize', 12);

subplot(1, 2, 2);
plot(y_mm, squeeze(T(cx, :, cz)), 'b-', 'LineWidth', 2.5);
grid on; xlabel('Lateral Distance y (mm)', 'FontSize', 11); ylabel('Temperature T (°C)', 'FontSize', 11);
title(sprintf('Lateral Profile T(y) at x = %.1f mm, z = %.2f mm', x_mm(cx), z_peak_mm), 'FontSize', 12);

sgtitle(sprintf('3D Thermal Lateral Diffusion Profiles - %s', source_label), ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(fig2, fullfile(run_folder, sprintf('%s_Lateral_Profiles_%s.png', source_label, timestamp_str)));

% --- Figure 3: Actual Diffused Temperature vs Wavelength ---
fig3 = figure('Name', 'Actual Diffused Temperature vs Wavelength', 'Color', 'w', 'Position', [100, 100, 900, 600]);
subplot(2, 1, 1);
plot(binned_wl, binned_pw * 1e3, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on; xlabel('Wavelength \lambda (nm)', 'FontSize', 11); ylabel('Optical Power (mW)', 'FontSize', 11);
title(sprintf('Input Optical Power Spectrum [%s]', source_label), 'FontSize', 12);

subplot(2, 1, 2);
plot(binned_wl, T_diffused_spectral, 'r-s', 'LineWidth', 2.2, 'MarkerFaceColor', 'r');
grid on; xlabel('Wavelength \lambda (nm)', 'FontSize', 11); ylabel('Diffused Temperature T (°C)', 'FontSize', 11);
title('Actual 3D PDE-Diffused Peak Temperature T_{diffused} vs Wavelength \lambda', 'FontSize', 12);

sgtitle(sprintf('Wavelength-Resolved Diffused Thermal Response: %s', source_label), ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(fig3, fullfile(run_folder, sprintf('%s_Spectral_Response_%s.png', source_label, timestamp_str)));

% --- Figure 4: FDM Convergence & Temperature Evolution Over Time & Iterations ---
fig4 = figure('Name', 'FDM Iterative Thermal Evolution & Convergence', 'Color', 'w', 'Position', [120, 120, 1200, 450]);

subplot(1, 3, 1);
plot(iter_history, temp_history, 'r-', 'LineWidth', 2.2);
grid on; xlabel('Iteration Number', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Peak Temperature T_{max} (°C)', 'FontSize', 11, 'FontWeight', 'bold');
title('(A) T_{max} vs. Iterations', 'FontSize', 12, 'FontWeight', 'bold');

subplot(1, 3, 2);
plot(time_history, temp_history, 'b-', 'LineWidth', 2.2);
grid on; xlabel('Elapsed Computation Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Peak Temperature T_{max} (°C)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('(B) T_{max} vs. Time (Final: %.2f s)', elapsed_sec), 'FontSize', 12, 'FontWeight', 'bold');

subplot(1, 3, 3);
semilogy(time_history, err_history, 'm-', 'LineWidth', 2.0);
grid on; xlabel('Elapsed Computation Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Pointwise Residual Error (°C)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('(C) Error vs. Time (Tol = 1e-6 °C in %d iters)', it), 'FontSize', 12, 'FontWeight', 'bold');

sgtitle(sprintf('FDM Thermal Evolution & Convergence Dynamics: %s', source_label), ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(fig4, fullfile(run_folder, sprintf('%s_FDM_Time_Convergence_%s.png', source_label, timestamp_str)));

% --- Figure 5: Optical Fluence Distribution (2D XZ Map & Axial Decay) ---
fig5 = figure('Name', 'Optical Fluence Distribution', 'Color', 'w', 'Position', [140, 140, 1100, 500]);

% Subplot 1: 2D XZ log10 Fluence Map
subplot(1, 2, 1);
fluence_xz_slice = squeeze(total_fluence(:, cy, :))'; % Depth z along rows, lateral x along cols
log_fluence_map = log10(max(fluence_xz_slice, 1e-12)); % Floor noise at -120 dB
max_log = max(log_fluence_map(:));
min_log = max_log - 4; % 4-decade dynamic range

imagesc(x_mm, z_mm, log_fluence_map);
axis image;
colormap(gca, parula);
cb = colorbar;
ylabel(cb, 'log_{10}(\Phi) [W/mm^2]', 'FontSize', 11);
clim([min_log, max_log]);
xlabel('Lateral Distance x (mm)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Tissue Depth z (mm)', 'FontSize', 11, 'FontWeight', 'bold');
title('(A) 2D XZ Fluence Map (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot 2: 1D Axial Depth Fluence Decay
subplot(1, 2, 2);
semilogy(z_mm, max(Phi_axial, 1e-12), 'k-', 'LineWidth', 2.3);
grid on;
xlabel('Tissue Depth z (mm)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Fluence Rate \Phi (W/mm^2)', 'FontSize', 11, 'FontWeight', 'bold');
title('(B) Axial Depth Fluence Decay \Phi(z)', 'FontSize', 12, 'FontWeight', 'bold');

sgtitle(sprintf('Optical Fluence Rate Distribution: %s', source_label), ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(fig5, fullfile(run_folder, sprintf('%s_Fluence_Map_%s.png', source_label, timestamp_str)));

fprintf('--> All output plots saved as PNG in:\n    %s\n\n', run_folder);

% =============================================================================
% HELPER FUNCTION: PIECEWISE TISSUE OPTICAL PROPERTIES (1/mm)
% =============================================================================
function [mu_a, mu_s, g, n] = get_tissue_optical_properties(wl_nm)
    n = 1.44; g = 0.95; 

    if (wl_nm >= 1400 && wl_nm <= 1500)
        mu_a = 2.80 * exp(-((wl_nm - 1450.0) / 22.0)^2);
    elseif (wl_nm >= 1850 && wl_nm <= 2030)
        mu_a = 12.0 * exp(-((wl_nm - 1940.0) / 30.0)^2);
    else
        mu_a = 0.0; % Zero absorption transparent windows
    end

    mu_s = 1.3295 * (wl_nm / 1000.0)^(-0.7665);
end

