clear variables;

% constants

tau = 0.1; % inversion time
ep = 1e-6; % tolerance parameter
Vm=5; % saturation voltage
Vc=1; %r emanent field
C=1;  % fitting parameter
a=1;  % fitting parameter
Cox = 1e-7;         % Oxide capacitance (F/cm^2)
phi = 0.7;

% Constants
mu = 1360;          % Electron mobility (cm^2/Vs)
W = 10e-8;          % Channel width (meters)
L = 1e-8;           % Channel length (meters)
Vth = 0.5;            % Threshold voltage (V)
lambda = 0.02;      % Channel-length modulation factor
Vds = 1;            % Fixed drain-source voltage (V)

V_G = linspace(0.1, 5,100);  % Gate voltage sweep (V)

%% iteration loop

for v_idx = 1:length(V_G)
    
    Vg = V_G(v_idx);

    % fprintf("V_gs = %d\n", Vg);
    V_mos = Vg;
    steps = 1;
    while true   
       
        V_fe = Vg - V_mos;

        Q_fe=(C/(2*a))*(atan((Vm+Vc)/a)-atan((Vm-Vc)/a))+(C/a)*(atan((V_fe-Vc)/a));
        Q_mos = Cox * (V_mos - phi);

        tolerence = abs(Q_fe-Q_mos);
        TOL(steps) = tolerence; 
        
        if tolerence < ep
            V_M(v_idx) = V_mos;
            fprintf("converged at %d\n for V_gs = %d\n", V_mos,Vg);            
            break;
        end

        % Update V_mos if tolerance error in large than epsilon
         V_mos = V_mos - 0.0001*tolerence;
         steps = steps + 1;
    end  
    
    

    % figure;
    % plot( V_M,TOL, 'b', 'LineWidth', 0.75);
    % xlabel('TOL');
    % ylabel('V_mos');
    % title('I-V Characteristics of a MOSFET (Ids vs Vgs)');
    % grid on;
end


%% Current vs V_gs graph -----------------

for i = 1:length(V_G)
    Vfe = V_M(i);
        if V_M(i) < Vth
            % Cut-off region, current is zero
            Ife(i) = 0;
        elseif V_M(i) >= Vth && Vds < (V_M(i) - Vth)
            % Linear region
            Ife(i) = mu * Cox * (W/L) * ((V_M(i)-Vth) * Vds - 0.5 * Vds^2)*(1+lambda*Vds) ;
        else
            % Saturation region
            Ife(i) = 0.5 * mu * Cox * (W/L) * (V_M(i)-Vth)^2 * (1 + lambda * Vds);
        end
    end

 for i = 1:length(V_G)
    if V_G(i) < Vth
        % Cut-off region, current is zero
        Ids(i) = 0;
    elseif V_G(i) >= Vth && Vds < (V_G(i) - Vth)
        % Linear region
        Ids(i) = mu * Cox * (W/L) * ((V_G(i) - Vth) * Vds - 0.5 * Vds^2)*(1+lambda*Vds) ;
    else
        % Saturation region
        Ids(i) = 0.5 * mu * Cox * (W/L) * (V_G(i) - Vth)^2 * (1 + lambda * Vds);
    end
end

% Plot the I-V Characteristics for ferroelectric fet
figure;
plot(V_G, Ids, 'b', 'LineWidth', 0.75);
hold on;

% Plot the I-V Characteristics

plot(V_G, Ife, 'red', 'LineWidth', 0.75);
xlabel('Vgs (V)');
ylabel('Drain Current Ids (A)');
title('I-V Characteristics of a MOSFET (Ids vs Vgs)');
grid on;

%% Current vs V_ds graph
V_g = 2;

Vds = 0:0.1:5;

Vm = 1.019528;

DIBL=1.127*(W\L);
% For FeFET--------------
for i = 1:length(Vds)
    if Vm < Vth
        % Cut-off region, current is zero
        Ife(i) = 0;
    elseif Vm >= Vth  && Vds(i) < (Vm - Vth)
        % Linear region
        Ife(i) = DIBL * ((Vm - Vth)* Vds(i) - 0.5 * Vds(i)^2)*(1+lambda*Vds(i));
    else
        % Saturation region
        Ife(i) = 0.5 * DIBL* (Vm- Vth)^2 * (1 + lambda * Vds(i));
    end
end

% For MOSFET-------------
for i = 1:length(Vds)
    if V_g < Vth
        % Cut-off region, current is zero
        Ids(i) = 0;
    elseif V_g >= Vth && Vds(i) < (V_g - Vth)
        % Linear region
        Ids(i) = DIBL * ((V_g - Vth)* Vds(i) - 0.5 * Vds(i)^2)*(1+lambda*Vds(i));
    else
        % Saturation region
        Ids(i) = 0.5 * DIBL* (V_g- Vth)^2 * (1 + lambda * Vds(i));
    end
end


% Plot the I-V Characteristics for ferroelectric fet
figure;
plot(Vds, Ids, 'b', 'LineWidth', 0.75);
hold on;

% Plot the I-V Characteristics

plot(Vds, Ife, 'red', 'LineWidth', 0.75);
xlabel('Vgs (V)');
ylabel('Drain Current Ids (A)');
title('I-V Characteristics of a MOSFET (Ids vs Vgs)');
grid on;

