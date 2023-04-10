

clear all; 
clc;

[t_inf,tb,p_inf,p_v,mul,rol,cc,lambdal,lw,alpha,rdot,ri,ro]=input();
% k=1e-6; %time-step
% N=10;
% n=N+2; % number of grid points 
% point = 1:1:n; 
% z = cos(0.5*pi - 0.5*pi*((point-1)/(n-1))); % collocation points 
% %Initial temeprature field and its first and second derivatives
% [T,Tz,Tzz] = initial(n,z,t_inf);

[h,t_r,t,r]=rk4_rpe(ri,rdot,mul,rol,p_v,p_inf,ro,cc);

plot(t/0.001,r(1,:)/ro,LineWidth=1.5,Color='black')

