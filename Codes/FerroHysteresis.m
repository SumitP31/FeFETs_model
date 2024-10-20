% A Ferroelectric Capacitor Mathematical Model for Spice Simulation
% CHAO-GANG WEI a , TIAN-LING REN a , JUN ZHU a &  LI-TIAN LIU a
clc 
clear variables
Vm=10;%saturation voltage
Vc=4;%remanent field
C=1; % fitting parameter
a=1; % fitting parameter
Vx=6; % reversal voltage for minor hysteresis 

% Primary Hysteresis Loop Equations for Q1, Q2, C1 and C2

V=-10:0.1:10;
for i=1:length(V)
    Q1(i)=(C/(2*a))*(atan((Vm+Vc)/a)-atan((Vm-Vc)/a))+(C/a)*(atan((V(i)-Vc)/a));
    Q2(i)=(C/(2*a))*(atan((Vm-Vc)/a)-atan((Vm+Vc)/a))+(C/a)*(atan((V(i)+Vc)/a));
    C1(i)=C/(a^2+(V(i)-Vc)^2);
    C2(i)=C/(a^2+(V(i)+Vc)^2);
end


% Minor Hysteresis Loop Equations for Q1, Q2, C1 and C2

Vn=-Vx:0.1:Vx;

% Equation for c'
Ct= C*((atan((Vm+Vc)/a)-atan((Vm-Vc)/a)+2*atan((Vx-Vc)/a))/(atan((Vx+Vc)/a)+atan((Vx-Vc)/a)));

for i=1:length(Vn)
    Qt1(i)=(Ct/(2*a))*(atan((Vx+Vc)/a)-atan((Vx-Vc)/a))+(Ct/a)*(atan((Vn(i)-Vc)/a));
    Qt2(i)=(Ct/(2*a))*(atan((Vx-Vc)/a)-atan((Vx+Vc)/a))+(Ct/a)*(atan((Vn(i)+Vc)/a));
    Cm1(i)=Ct/(a^2+(V(i)-Vc)^2);
    Cm2(i)=Ct/(a^2+(V(i)+Vc)^2);
end



% Primary Hysteresis Loop Plot
figure
plot(V,Q1,'b','Marker','square')
hold on
plot(V,Q2,'r','Marker','diamond')
plot(V,C1,'g','LineWidth',2)
plot(V,C2,'y','LineWidth',2)

title("Primary Hysteresis Loop ");
xlabel("Voltage (Volts)")
ylabel("Q,C")
legend("Primary Lower Branch","Primary Upper Branch","Capacitance Lower Branch","Capacitance Upper Branch")

% Minor Hysteresis Loop Plot
figure
plot(V,Q1,'b','Marker','square') % from primary loop eq
hold on
plot(V,Q2,'r','Marker','diamond') % from primary loop eq
plot(Vn,Qt1,'g','LineWidth',2)
plot(Vn,Qt2,'black','LineWidth',2)

title("Minor Hysteresis Loop ");
xlabel("Voltage (Volts)")
ylabel("Q,C")
legend("Primary Lower Branch","Primary Upper Branch","Minor Lower Branch","Minor Upper Branch")

