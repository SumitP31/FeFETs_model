% MATLAB Code for FeFET Model using Preisach Theory

% Constants
VDD = 1;       % Supply Voltage (V)
tau = 1;    % Time constant for dynamic switching behavior
error_tolerance = 1e-1;   % Convergence criterion

% BSIM Model Parameters (for MOSFET part)
beta = 0.02;   % Example value for BSIM parameter (V/A)

% Ferroelectric Layer Parameters
alpha_center = -1; % Coercive field lower bound
beta_center = 1;   % Coercive field upper bound
N_domains = 10;  % Number of domains in the Preisach model
epsilon_Fe = 1e-9; % Dielectric constant of ferroelectric layer (F/m)

% Generate Random Coercive Fields for Each Domain (Normal Distribution)
alpha_i = normrnd(alpha_center, 0.1, [N_domains, 1]);
beta_i = normrnd(beta_center, 0.1, [N_domains, 1]);

% Initialize Variables
V_G = linspace(-2, 2, 10);  % Gate voltage sweep (V)
I_D = zeros(size(V_G));      % Drain current (A)

% Preisach Model for FE Layer and Iterative Solution
for v_idx = 1:length(V_G)
    Vg = V_G(v_idx);
    Vmos = 0; % Initial guess for V_MOS
    
    % Iterative solution for each gate voltage
    while true
        Qmos = beta * Vmos;  % Using BSIM model for MOSFET (simplified)
        Vfe = Vg - Vmos;     % Voltage drop across FE layer

        % Calculate charge in FE layer using Preisach model
        P_i = zeros(N_domains, 1);  % Polarization state of each domain
        E = Vfe / epsilon_Fe;       % Electric field in FE layer
        
        for i = 1:N_domains
            if E <= alpha_i(i)
                P_i(i) = -1;   % Polarization down state
            elseif E >= beta_i(i)
                P_i(i) = 1;    % Polarization up state
            else
                P_i(i) = 0;    % No change in polarization
            end
        end

        % Total spontaneous polarization of FE layer
        Qfe = sum(P_i) / N_domains;

        % Dynamic behavior of FE layer (RC-like model)
        Veff = Vfe * exp(-tau * (Vfe - Qfe));

        % Check convergence
        if abs(Qfe - Qmos) < error_tolerance
            break;
        end
        
        % Update Vmos for next iteration
        Vmos = Vmos + 0.1 * (Qfe - Qmos);
    end
    
    % Calculate drain current based on the final Vmos
    I_D(v_idx) = Qmos * (Vg - Vmos);
end

 % Plot the Gate Voltage vs Drain Current
 figure;
 plot(V_G, I_D, 'LineWidth', 2);
 xlabel('Gate Voltage (V)');
 ylabel('Drain Current (A)');
 title('FeFET Gate Voltage vs Drain Current');
 grid on;
