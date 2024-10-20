clear all;
hfin=1;
tsi=20
W=2*hfin*tsi;
Vth=1;
Leff=180;
DIBL=1.127*(W\Leff);
vth=Vth-DIBL;
lambda=25*10^-5;
u=1360;
ld=0.1;
Ec=1.5*106;
Cox=60;
Vds=0:0.5:10;
vgs=input('ENTER THE Vgs in volts');
m=length(Vds);

current = zeros(1, m);
current1 = zeros(1, m);

for i=1:m

if vgs < Vth
current(1,i)=0;
current1(1,i)=0;
elseif Vds(i) >= (vgs - Vth)
current(1,i)=0.5* DIBL * ((vgs - Vth)^2)*(1+lambda*Vds(i));%Added DIBL
elseif Vds(i) < (vgs - Vth)
current(1,i)= DIBL*((vgs-Vth)*Vds(i) - 0.5*(Vds(i)^2))*(1+lambda*Vds(i)); %Simplified equation by approximation
end

end
plot(Vds,current(1,:),'b')
xlabel('Vds, V')
ylabel('Drain Current,A')
title('I-V Characteristics of a MOSFET')