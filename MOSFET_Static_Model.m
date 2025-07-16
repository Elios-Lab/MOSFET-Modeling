clear; clc;

% --- Realistic parameters for IRF740 ---
Vth = 4.7;               % Typical threshold voltage [V]
KPLIN = 2.4;             % Gain in the linear region [A/V^2]
KPSAT = 8.84;            % Gain in the saturation region [A/V^2]
theta = 0.3;             % Early effect coefficient

% --- Sweep VDS and VGS ---
VDS = linspace(0, 100, 1000);             % Up to 100 V
VGS_vals = -20:1:20;                      % Full range for IRF740
colors = lines(length(VGS_vals));

figure;

%% --- SUBPLOT 1: IDS vs VDS (linear + saturation regions) ---
subplot(2,1,1); hold on; grid on;
title('IRF740 - I_{DS} vs V_{DS} for different V_{GS}');
xlabel('V_{DS} [V]'); ylabel('I_{DS} [A]');

for idx = 1:length(VGS_vals)
    VGS = VGS_vals(idx);
    if VGS <= Vth
        continue; % Cutoff region, skip plotting
    end
    IDS = zeros(size(VDS));
    VDS_border = VGS - Vth;
    
    % Current at the boundary (for continuity with saturation region)
    ID_border = (KPLIN * ((VGS - Vth) * VDS_border - (0.5 * KPLIN / KPSAT * VDS_border^2))) / ...
                (1 + theta * VDS_border);

    for i = 1:length(VDS)
        vds = VDS(i);
        if vds < VDS_border
            IDS(i) = (KPLIN * ((VGS - Vth) * vds - (0.5 * KPLIN / KPSAT * vds^2))) / ...
                     (1 + theta * VDS_border);  % Linear region
        else
            IDS(i) = ID_border * (1 + theta * (vds - VDS_border));  % Continuous saturation
        end
    end

    IDS = min(IDS, 10);  % Limit to 10 A (from datasheet)
    
    plot(VDS, IDS, 'LineWidth', 2, 'Color', colors(mod(idx-1, size(colors,1))+1, :), ...
        'DisplayName', sprintf('V_{GS} = %.1f V', VGS));
end
legend show;
ylim([0 12]); xlim([0 100]);

%% --- SUBPLOT 2: IDS vs VGS (saturation region only) ---
subplot(2,1,2); hold on; grid on;
title('Saturation region: Parabolic behavior of I_{DS} vs V_{GS}');
xlabel('V_{GS} [V]'); ylabel('I_{DS} [A]');

VGS_range = linspace(-20, 20, 500);
VDS_fixed = 100;  % Fixed to ensure saturation condition

IDS_sat = zeros(size(VGS_range));
for i = 1:length(VGS_range)
    VGS = VGS_range(i);
    if VGS <= Vth
        IDS_sat(i) = 0;
    else
        IDS_sat(i) = (KPSAT / 2 * (VGS - Vth)^2) / (1 + theta * VDS_fixed);
    end
end

IDS_sat = min(IDS_sat, 10);  % Datasheet limit
plot(VGS_range, IDS_sat, 'b', 'LineWidth', 2);
legend(sprintf('V_{DS} = %.1f V (saturation)', VDS_fixed));
ylim([0 12]); xlim([-5 20]);  % Zoom into useful region

%
% The two formulas are not constructed to perfectly match at the 
% transition point VDS = VGS − Vth.
%
% Specifically:
%
% A) The model is not "continuous by design": the end of the linear region 
% does not numerically match the start of the saturation region.
% B) The two coefficients KPLIN and KPSAT are different, so even though 
% the terms look similar, they yield different values at the transition 
% point → hence the visible offset in the plot.
%
% To observe continuity, you must force KPLIN = KPSAT, adjusting 
% the current formula accordingly.
