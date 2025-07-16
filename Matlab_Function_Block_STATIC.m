% === Model Parameters ===
Vth   = 4.7;        % Threshold voltage [V]
KPLIN = 2.4;        % Gain in the linear (triode) region [A/V^2]
KPSAT = 11.7;       % Gain in the saturation region [A/V^2]
theta = 0.3;        % Early effect coefficient [1/V]

% === IDS Current Calculation ===
if vgs <= Vth
    % Cut-off region (transistor is off)
    ids = 0;

elseif vds < (vgs - Vth)
    % Linear (triode) region
    ids = KPLIN * ((vgs - Vth) * vds - 0.5 * vds^2) / (1 + theta * vds);

else
    % Saturation region
    ids = 0.5 * KPSAT * (vgs - Vth)^2 / (1 + theta * vds);
end
