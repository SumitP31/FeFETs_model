% MOSFET I-V Characteristics (Ids vs Vgs) - MATLAB Code

clear all;
clc;

% Constants
mu = 1360;          % Electron mobility (cm^2/Vs)
W = 10e-8;          % Channel width (meters)
L = 1e-8;           % Channel length (meters)
Cox = 1e-7;         % Oxide capacitance (F/cm^2)
Vth = 0.5;            % Threshold voltage (V)
lambda = 0.02;      % Channel-length modulation factor
Vds = 1;            % Fixed drain-source voltage (V)

% Vgs range for the plot
Vgs = 0:0.1:2;      % Gate-source voltage range (V)

% Preallocate current array
Ids = zeros(1, length(Vgs));

% Loop through the Vgs values to compute Ids
for i = 1:length(Vgs)
    if Vgs(i) < Vth
        % Cut-off region, current is zero
        Ids(i) = 0;
    elseif Vgs(i) >= Vth && Vds < (Vgs(i) - Vth)
        % Linear region
        Ids(i) = mu * Cox * (W/L) * ((Vgs(i) - Vth) * Vds - 0.5 * Vds^2) ;
    else
        % Saturation region
        Ids(i) = 0.5 * mu * Cox * (W/L) * (Vgs(i) - Vth)^2 * (1 + lambda * Vds);
    end
end

% Plot the I-V Characteristics
figure;
plot(Vgs, Ids, 'b', 'LineWidth', 0.75);
xlabel('Vgs (V)');
ylabel('Drain Current Ids (A)');
title('I-V Characteristics of a MOSFET (Ids vs Vgs)');
grid on;
