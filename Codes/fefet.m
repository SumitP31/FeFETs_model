% FeFET Characteristics - MATLAB Code

clear all;
clc;

% Constants
mu = 1360;          % Mobility (cm^2/Vs)
W = 10e-6;          % Channel width (meters)
L = 1e-6;           % Channel length (meters)
Cox = 1e-7;         % Oxide capacitance (F/cm^2)
Cfe = 2e-7;         % Ferroelectric capacitance (F/cm^2)
Vth0 = 1;           % Intrinsic threshold voltage (V)
Pr = 30e-6;         % Remnant polarization (C/cm^2)
lambda = 0.02;      % Channel-length modulation factor
Vds = 0:0.1:5;      % Drain-source voltage range (V)

% Gate Voltage Input
Vgs_ext = input('Enter the external gate voltage Vgs (in volts): ');

% Calculate effective gate capacitance
Cg = (Cfe * Cox) / (Cfe + Cox);

% Calculate polarization and effective threshold voltage
Vth = Vth0 - (Pr / Cox);

% Preallocate current array
Ids_linear = zeros(1, length(Vds));
Ids_saturation = zeros(1, length(Vds));

% Loop through the Vds values to compute Ids
for i = 1:length(Vds)
    if Vds(i) < (Vgs_ext - Vth)
        % Linear region
        Ids_linear(i) = mu * Cg * (W/L) * ((Vgs_ext - Vth) * Vds(i) - 0.5 * Vds(i)^2);
    else
        % Saturation region
        Ids_saturation(i) = 0.5 * mu * Cg * (W/L) * (Vgs_ext - Vth)^2 * (1 + lambda * Vds(i));
    end
end

% Plot the I-V Characteristics
figure;
plot(Vds, Ids_linear, 'b', Vds, Ids_saturation, 'r--');
xlabel('Vds (V)');
ylabel('Drain Current Ids (A)');
title('I-V Characteristics of FeFET');
legend('Linear Region', 'Saturation Region');
grid on;
