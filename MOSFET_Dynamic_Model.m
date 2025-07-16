clear; clc;

% --- Static model parameters ---
Vth = 4.7;               % Threshold voltage [V]
KPLIN = 2.4;             % Linear region gain [A/V^2]
KPSAT = 8.84;            % Saturation region gain [A/V^2]
theta = 0.3;             % Early effect

% --- Parasitic capacitive parameters ---
Cgs = 1373e-12;          % [F]
Cgd = 27e-12;            % [F]
Cds = 193e-12;           % [F] - Does not directly affect IDS

% --- Current limits ---
Imax_continuous = 10;    % [A]
Imax_pulse = 40;         % [A]
t_max_pulse = 100e-6;    % [s]

% --- Time and signals ---
t_max = 10e-6;                           % 10 µs
t = linspace(0, t_max, 5000);
dt = t(2) - t(1);

f_vgs = 1e6;                             % 1 MHz
f_vds = 0.5e6;                           % 0.5 MHz

VGS_t = 10 * sin(2*pi*f_vgs*t);          % VGS: 10 Vpp @ 1 MHz
VDS_t = 200 + 100 * sin(2*pi*f_vds*t + pi/4); % VDS: 200V offset + 100 Vpp @ 0.5 MHz

dVGS_dt = gradient(VGS_t, dt);
dVDS_dt = gradient(VDS_t, dt);

% --- Current ---
IDS_static  = zeros(size(t));
IDS_dynamic = zeros(size(t));

for i = 1:length(t)
    VGS = VGS_t(i);
    VDS = VDS_t(i);

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

    % Static clipping
    IDS_static(i) = min(IDS_static(i), Imax_continuous);

    % Dynamic current
    IDS_dynamic(i) = IDS_static(i) + ...
                     (Cgs + Cgd) * dVGS_dt(i) - Cgd * dVDS_dt(i);
end

% --- Dynamic limits ---
pulse_in_progress = false;
pulse_duration = 0;

for i = 1:length(t)
    if IDS_dynamic(i) > Imax_continuous
        if ~pulse_in_progress
            pulse_in_progress = true;
            pulse_duration = dt;
        else
            pulse_duration = pulse_duration + dt;
        end

        if IDS_dynamic(i) > Imax_pulse || pulse_duration > t_max_pulse
            IDS_dynamic(i) = Imax_continuous;
        end
    else
        pulse_in_progress = false;
        pulse_duration = 0;
    end
end

%% --- Plot 1: I_DS(t) ---
figure; hold on; grid on;
title('IRF740 - I_{DS}(t) Current with Dynamic Effects and Limits');
plot(t*1e6, IDS_static, 'b--', 'DisplayName', 'Static I_{DS}');
plot(t*1e6, IDS_dynamic, 'r', 'LineWidth', 2, 'DisplayName', 'Dynamic I_{DS}');
yline(Imax_continuous, 'k--', 'LineWidth', 1.5, 'DisplayName', '10 A Limit');
yline(Imax_pulse, 'm--', 'LineWidth', 1.2, 'DisplayName', '40 A Limit');
xlabel('Time [\mus]');
ylabel('I_{DS} [A]');
legend;

%% --- Plot 2: I_DS vs V_DS ---
figure; hold on; grid on;
title('IRF740 - I_{DS}(V_{DS}) Current with Dynamics');
plot(VDS_t, IDS_dynamic, 'r', 'LineWidth', 1.5, 'DisplayName', 'Dynamic I_{DS}');
plot(VDS_t, IDS_static, 'b--', 'LineWidth', 1, 'DisplayName', 'Static I_{DS}');
yline(Imax_continuous, 'k--', 'LineWidth', 1.5, 'DisplayName', '10 A Limit');
yline(Imax_pulse, 'm--', 'LineWidth', 1.2, 'DisplayName', '40 A Limit');
xlabel('V_{DS} [V]');
ylabel('I_{DS} [A]');
legend;

%% --- Plot 3: V_GS and V_DS vs Time with Conduction Zones Highlighted ---
figure; hold on; grid on;
title('V_{GS} and V_{DS} with Highlighted Conduction Zones');
xlabel('Time [\mus]');
ylabel('Voltage [V]');
plot(t*1e6, VGS_t, 'b', 'LineWidth', 1.5, 'DisplayName', 'V_{GS}');
plot(t*1e6, VDS_t, 'r', 'LineWidth', 1.5, 'DisplayName', 'V_{DS}');
yline(Vth, 'k--', 'DisplayName', 'V_{th}');

% Highlight conduction zones (VGS > Vth)
hold on;
idx_on = VGS_t > Vth;
on_regions = bwconncomp(idx_on);  % Identify ON zones

for k = 1:on_regions.NumObjects
    idxs = on_regions.PixelIdxList{k};
    t_fill = t(idxs) * 1e6;
    fill([t_fill fliplr(t_fill)], ...
         [min(VDS_t)*ones(size(t_fill)), max(VDS_t)*ones(size(t_fill))], ...
         [0.3 1 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.6);   % Light green
end

legend;
