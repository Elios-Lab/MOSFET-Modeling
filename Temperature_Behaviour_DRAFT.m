clear; clc;

% --- Static Model Parameters ---
Vth_25 = 4.7;             % Threshold voltage at 25°C [V]
KPLIN = 2.4;              % Linear region gain [A/V^2]
KPSAT = 8.84;             % Saturation region gain [A/V^2]
theta = 0.3;              % Early effect coefficient

% --- Parasitic Capacitances ---
Cgs = 1373e-12;           % Gate-source capacitance [F]
Cgd = 27e-12;             % Gate-drain capacitance [F]
Cds = 193e-12;            % Drain-source capacitance [F] – doesn't directly affect IDS

% --- Thermal Parameters ---
Tamb = 25;                % Ambient temperature [°C]
Rth_ja = 62;              % Thermal resistance junction-to-ambient [°C/W]
Cth = 0.5;                % Thermal capacitance [J/°C]
Rds_25 = 0.48;            % Rds(on) at 25°C [Ohm]
alpha = 0.005;            % Rds(on) temp coefficient [1/°C]
beta = 2.5e-3;            % Vth temperature coefficient [V/°C]
Tj_max = 150;             % Max junction temperature [°C]

% --- Current Limits ---
Imax_continuous = 10;     % Max continuous current [A]
Imax_pulse = 40;          % Max pulse current [A]
t_pulse_max = 100e-6;     % Max pulse duration [s]

% --- Time and Signals ---
t_max = 1;                            % Total simulation time [s]
n_points = 1e5;                       % Number of time steps
t = linspace(0, t_max, n_points);
dt = t(2) - t(1);

f_vgs = 1e3;                          % VGS frequency: 1 kHz
f_vds = 0.5e3;                        % VDS frequency: 0.5 kHz

VGS_t = 10 * sin(2*pi*f_vgs*t);       % VGS: 10 Vpp @ 1 kHz
VDS_t = 200 + 100 * sin(2*pi*f_vds*t + pi/4); % VDS: 200V offset + 100 Vpp @ 0.5 kHz

dVGS_dt = gradient(VGS_t, dt);
dVDS_dt = gradient(VDS_t, dt);

% --- Current and Temperature ---
IDS_static  = zeros(size(t));
IDS_dynamic = zeros(size(t));
Tj_t        = zeros(size(t));           % Junction temperature
Rds_t       = zeros(size(t));
Vth_t       = zeros(size(t));

% Initialize temperature
Tj_t(1) = Tamb;

for i = 1:length(t)
    VGS = VGS_t(i);
    VDS = VDS_t(i);

    % --- Dynamic Junction Temperature (RC Model) ---
    if i == 1
        Tj = Tamb;
    else
        P_diss = IDS_dynamic(i-1)^2 * Rds_t(i-1);  % Power loss
        dT = dt / Cth * (P_diss - (Tj_t(i-1) - Tamb) / Rth_ja);
        Tj = Tj_t(i-1) + dT;
    end
    Tj_t(i) = Tj;

    % Thermal shutdown protection
    if Tj > Tj_max
        IDS_static(i) = 0;
        IDS_dynamic(i) = 0;
        Rds_t(i) = Rds_25 * (1 + alpha * (Tj - 25));
        Vth_t(i) = Vth_25 - beta * (Tj - 25);
        continue;  % Skip to next time step
    end

    % --- Update temperature-dependent parameters ---
    Rds = Rds_25 * (1 + alpha * (Tj - 25));
    Vth = Vth_25 - beta * (Tj - 25);
    Rds_t(i) = Rds;
    Vth_t(i) = Vth;

    % --- Static Current Calculation ---
    if VGS <= Vth
        IDS_static(i) = 0;
    elseif VDS < VGS - Vth
        IDS_static(i) = (KPLIN * ((VGS - Vth) * VDS - ...
                         (0.5 * KPLIN / KPSAT * VDS^2))) / ...
                         (1 + theta * VDS);
    else
        VDS_border = VGS - Vth;
        ID_border = (KPLIN * ((VGS - Vth) * VDS_border - ...
                     (0.5 * KPLIN / KPSAT * VDS_border^2))) / ...
                     (1 + theta * VDS_border);
        IDS_static(i) = ID_border * (1 + theta * (VDS - VDS_border));
    end

    % Clamp static current
    IDS_static(i) = min(IDS_static(i), Imax_continuous);

    % --- Dynamic Current (includes capacitive effects) ---
    IDS_dynamic(i) = IDS_static(i) + ...
                     (Cgs + Cgd) * dVGS_dt(i) - Cgd * dVDS_dt(i);
end

% --- Pulse current limiting ---
in_pulse = false;
pulse_duration = 0;

for i = 1:length(t)
    if IDS_dynamic(i) > Imax_continuous
        if ~in_pulse
            in_pulse = true;
            pulse_duration = dt;
        else
            pulse_duration = pulse_duration + dt;
        end

        if IDS_dynamic(i) > Imax_pulse || pulse_duration > t_pulse_max
            IDS_dynamic(i) = Imax_continuous;
        end
    else
        in_pulse = false;
        pulse_duration = 0;
    end
end

%% --- Plot 1: IDS vs Time ---
figure; hold on; grid on;
title('IRF740 - I_{DS}(t) with dynamic effects and limits');
plot(t, IDS_static, 'b--', 'DisplayName', 'Static I_{DS}');
plot(t, IDS_dynamic, 'r', 'LineWidth', 2, 'DisplayName', 'Dynamic I_{DS}');
yline(Imax_continuous, 'k--', 'LineWidth', 1.5, 'DisplayName', '10 A Limit');
yline(Imax_pulse, 'm--', 'LineWidth', 1.2, 'DisplayName', '40 A Limit');
xlabel('Time [s]');
ylabel('I_{DS} [A]');
legend;

%% --- Plot 2: IDS vs VDS ---
figure; hold on; grid on;
title('IRF740 - I_{DS}(V_{DS}) with dynamic effects');
plot(VDS_t, IDS_dynamic, 'r', 'LineWidth', 1.5, 'DisplayName', 'Dynamic I_{DS}');
plot(VDS_t, IDS_static, 'b--', 'LineWidth', 1, 'DisplayName', 'Static I_{DS}');
yline(Imax_continuous, 'k--', 'LineWidth', 1.5, 'DisplayName', '10 A Limit');
yline(Imax_pulse, 'm--', 'LineWidth', 1.2, 'DisplayName', '40 A Limit');
xlabel('V_{DS} [V]');
ylabel('I_{DS} [A]');
legend;

%% --- Plot 3: Junction Temperature over Time ---
figure; hold on; grid on;
title('Junction Temperature T_j(t)');
plot(t, Tj_t, 'm', 'LineWidth', 1.5);
yline(Tj_max, 'k--', 'DisplayName', 'T_{j max} = 150°C');
xlabel('Time [s]');
ylabel('T_j [°C]');
legend;
