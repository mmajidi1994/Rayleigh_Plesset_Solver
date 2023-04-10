function [t_inf,tb,p_inf,p_v,mul,rol,cc,lambdal,lw,alpha,rdot,ri,ro]=input() 
%%%%%
t_inf=300; %ambient temperature
tb=t_inf;
p_inf=100*1000; %ambient pressure
p_v=4239.651;  %vapour pressure 
mul=0.000797;  %dynamic viscosity of water 
rol=995.65;    %water density
cc=0.0712;     %surface tension 

%%%%%%
lambdal=0.6; %water heat conductivity 
lw=334000; %latent heat of water
alpha= 0.6; %thermal diffusivity of water
%%%%%%
rdot=0; %Initial velocity of the bubble interface
ri=0.0004; %initial radius 
ro=0.0002; %equilibrium radius
%%%%%

end 
