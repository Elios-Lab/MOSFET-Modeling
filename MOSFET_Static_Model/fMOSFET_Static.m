
function ids = fMOSFET_Static(vgs, vds)
% Static MOSFET model based on article
Vth = 4.7;
KPLIN = 2.4;
KPSAT = 8.84;
theta = 0.3;

if vgs <= Vth
    ids = 0; % Cut-off
elseif vds < (vgs - Vth)
    ids = KPLIN * ((vgs - Vth) * vds - 0.5 * vds^2); % Linear region
else
    ids = KPSAT * (vgs - Vth)^2 * (1 + theta * vds); % Saturation
end
end
