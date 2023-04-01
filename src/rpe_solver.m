
clc
clear all

 
ri=0.0004; %initial radius 
ro=0.0002; %equilibrium radius


[t,y] = ode45(@vdp1,[0 0.0005],[ri; 0]);
plot(t/0.001,y(:,1)/ro,'-o')
title('Solution of simpleset Rayleigh Plesset equation (when bubble temperature is constant)');
xlabel('Time t');
ylabel('Solution y');
legend('R')
grid on

function dydt = vdp1(t,y)

p_inf=100*1000; %100kpa
p_v=4239.651;

mul=0.000797;
rol=995.65;
cc=0.0712;

ri=0.0004; %initial radius 
ro=0.0002; %equilibrium radius

tem= -1.5*((y(2).^2)/y(1))-4*(mul/rol)*(y(2)/(y(1).^2))...
    +(1/(rol*y(1)))*(p_v-p_inf+p_inf*(ro/y(1)).^3 -2*cc/y(1));

dydt = [y(2); tem ];

end 
