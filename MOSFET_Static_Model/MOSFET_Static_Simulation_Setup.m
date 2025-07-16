
% Script per la simulazione del MOSFET statico in Simulink
% Deve essere accompagnato da un modello Simulink che include:
% - Il blocco MATLAB Function con fMOSFET_Static
% - Controlled Current Source pilotato dalla corrente calcolata
% - RL Load, sorgente VDS, generatore di impulsi VGS

% Parametri di simulazione
VDD = 150;         % Tensione Drain-Source
Rload = 69;        % Resistenza del carico
Lload = 0.7e-6;    % Induttanza del carico
RG = 15;           % Resistenza di gate
VGS_on = 10;       % Impulso VGS ON
VGS_off = 0;       % VGS OFF

% Pulse Generator per VGS (Simulink)
pulse_amp = VGS_on;
pulse_period = 2e-6;  % 2 us periodo
pulse_width = 50;     % % duty cycle
pulse_delay = 0;

% Tempo di simulazione
t_end = 5e-6; % 5 microsecondi
