function IDS = mosfet_ids_model(VGS, VDS, dVGS_dt, dVDS_dt)
%#codegen
% Computes IDS based on VGS, VDS, and their derivatives, using the IRF740 dynamic model

% --- IRF740 model parameters ---
Vth = 4.7;               % Threshold voltage [V]
KPLIN = 2.4;             % Gain in the linear region [A/V^2]
KPSAT = 8.84;            % Gain in the saturation region [A/V^2]
theta = 0.3;             % Early effect coefficient

% --- Capacitive parameters ---
Cgs = 1373e-12;          % Gate-to-source capacitance [F]
Cgd = 27e-12;            % Gate-to-drain capacitance [F]

% --- Current limits ---
Imax_continuous = 10;    % Maximum continuous current [A]
Imax_pulse      = 40;    % Maximum pulse current [A]
% Note: pulse duration should be handled externally in Simulink if needed

% --- Static current ---
if VGS <= Vth
    IDS_static = 0;
elseif VDS < VGS - Vth
    IDS_static = (KPLIN * ((VGS - Vth) * VDS - ...
                 (0.5 * KPLIN / KPSAT * VDS^2))) / ...
                 (1 + theta * VDS);
else
    VDS_border = VGS - Vth;
    ID_border = (KPLIN * ((VGS - Vth) * VDS_border - ...
                 (0.5 * KPLIN / KPSAT * VDS_border^2))) / ...
                 (1 + theta * VDS_border);
    IDS_static = ID_border * (1 + theta * (VDS - VDS_border));
end

% --- Static clipping ---
IDS_static = min(IDS_static, Imax_continuous);

% --- Dynamic terms (capacitive components) ---
IDS_dynamic = IDS_static + (Cgs + Cgd) * dVGS_dt - Cgd * dVDS_dt;

% --- Pulse clipping (duration must be managed externally) ---
IDS = min(IDS_dynamic, Imax_pulse);
end
