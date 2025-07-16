
% Script per la simulazione del MOSFET con modello dinamico (capacità inter-elettrodo)
% Aggiunge CGS, CGD e CDS tra i terminali G, D e S come indicato nell'articolo

% Capacità inter-elettrodo dal paper (IRF740)
CGS = 1373e-12;  % Gate-Source capacitance [F]
CGD = 27e-12;    % Gate-Drain capacitance [F]
CDS = 193e-12;   % Drain-Source capacitance [F]

% Inserire questi valori nei blocchi "Capacitor" di Simscape
% - CGS: collegare tra Gate e Source
% - CGD: collegare tra Gate e Drain
% - CDS: collegare tra Drain e Source

% Il resto del circuito (Pulse Generator, RL load, fMOSFET_Static) rimane invariato
% Simulare con stesso setup di tempi: t_end = 5e-6
