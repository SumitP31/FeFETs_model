% A Ferroelectric Capacitor Mathematical Model for Spice Simulation
% CHAO-GANG WEI a , TIAN-LING REN a , JUN ZHU a &  LI-TIAN LIU a
clc 
clear all
Vm=8;%saturation voltage
Vc=3;%remanent field
C=2;
a=0.9;
x1=-8:0.1:8;
x2=7.9:-0.1:-6;
x3=-5.9:0.1:4;
x4=3.9:-0.1:-3;
x5=-2.9:0.1:2;
x6=1.9:-0.01:-1;
V=[x1 x2 x3 x4 x5];
%x=-pi:0.1:0;
%V=3*sin(x);
%x=-5:0.1:5
y=gradient(V);

% Vx=5;
% Vn=-Vx:0.1:Vx;
Vn=V(length(x1));
Vx=Vn;
[Q1,Q2,C1,C2]=satu(Vn,a,Vc,Vm,C);
Vc2=Vc;
Vc1=Vc;
Ct=C*((atan((Vm+Vc)/a)-atan((Vm-Vc)/a)+2*atan((Vx-Vc)/a))/(atan((Vx+Vc)/a)+atan((Vx-Vc)/a)));
for i=1:length(V)
    
    if  i<length(V) && (V(i+1)-V(i))>0
        Vc=Vc1;
        Qt1(i)=(Ct/(2*a))*(atan((Vx+Vc)/a)-atan((Vx-Vc)/a))+(Ct/a)*(atan((V(i)-Vc)/a));
        if i>1 && i<length(V) && Qt1(i)*Qt1(i-1)*100<0
            Vc2=abs(V(i));
            % C=Ct;
            %  Ct=C*((atan((Vm+Vc)/a)-atan((Vm-Vc)/a)+2*atan((Vx-Vc)/a))/(atan((Vx+Vc)/a)+atan((Vx-Vc)/a)));
        end
    end
    if  i<length(V) && (V(i+1)-V(i))<0
        Vc=Vc2;
        Qt1(i)=(Ct/(2*a))*(atan((Vx-Vc)/a)-atan((Vx+Vc)/a))+(Ct/a)*(atan((V(i)+Vc)/a));
        if i>1 && i<length(V) && Qt1(i)*Qt1(i-1)*100<0
            Vc1=abs(V(i));
             %  C=Ct;
             % Ct=C*((atan((Vm+Vc1)/a)-atan((Vm-Vc1)/a)+2*atan((Vx-Vc1)/a))/(atan((Vx+Vc1)/a)+atan((Vx-Vc1)/a)));
        end
    end
    if  i==length(V) && (V(i)-V(i-1))>0
        Qt1(i)=(Ct/(2*a))*(atan((Vx+Vc)/a)-atan((Vx-Vc)/a))+(Ct/a)*(atan((V(i)-Vc)/a));
    end
    if  i==length(V) && (V(i)-V(i-1))<0
        Qt1(i)=(Ct/(2*a))*(atan((Vx-Vc2)/a)-atan((Vx+Vc2)/a))+(Ct/a)*(atan((V(i)+Vc2)/a));
    end
    if  i<length(V)-1 && (V(i+1)-V(i))*(V(i+2)-V(i+1))<0
        Vm=Vx;
        Vx=abs(V(i+1));
        C=Ct;
        Ct=C*((atan((Vm+Vc)/a)-atan((Vm-Vc)/a)+2*atan((Vx-Vc)/a))/(atan((Vx+Vc)/a)+atan((Vx-Vc)/a)));

    end
end
figure
plot(V,Qt1,'black','LineWidth',0.5,'Marker','.')
hold on
title("Hysteresis Loop");
xlabel("Voltage (Volts)")
ylabel("Q")
figure
plot(V,'blue','LineWidth',0.5,'Marker','.')
title("Applied Voltage");
%xlabel("Step size")
ylabel("Voltage (Volts)")

function [Q1,Q2,C1,C2]=satu(Vn,a,Vc,Vm,C)
    for i=1:length(Vn)
        Q1(i)=(C/(2*a))*(atan((Vm+Vc)/a)-atan((Vm-Vc)/a))+(C/a)*(atan((Vn(i)-Vc)/a));
        Q2(i)=(C/(2*a))*(atan((Vm-Vc)/a)-atan((Vm+Vc)/a))+(C/a)*(atan((Vn(i)+Vc)/a));
        C1(i)=C/(a^2+(Vn(i)-Vc)^2);
        C2(i)=C/(a^2+(Vn(i)+Vc)^2);
    end
end
