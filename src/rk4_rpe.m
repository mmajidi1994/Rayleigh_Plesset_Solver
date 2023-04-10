function [h,t_r,t,r]=rk4_rpe(ri,rdot,mul,rol,p_v,p_inf,ro,cc)

h=1e-8;               % step size
t_r=0.0005;
t = 0:h:t_r;         
r = zeros(2,length(t)); 
% initial condition
r(1,1) =ri;  
r(2,1) =rdot;  

F_tr1 = @(t,r1,r2) r2;               % change the function as you desire
F_tr2 = @(t,r1,r2) -1.5*((r2.^2)/r1)-4*(mul/rol)*(r2/(r1.^2))...
     +(1/(rol*r1))*(p_v-p_inf+p_inf*(ro/r1).^3 -2*cc/r1); % change the function as you desire

for i=1:(length(t)-1)                              % calculation loop
    k_11 = F_tr1(t(i),r(1,i),r(2,i));
    k_12 = F_tr2(t(i),r(1,i),r(2,i));
    k_21 = F_tr1(t(i)+0.5*h,r(1,i)+0.5*h*k_11,r(2,i)+0.5*h*k_12);
    k_22 = F_tr2(t(i)+0.5*h,r(1,i)+0.5*h*k_11,r(2,i)+0.5*h*k_12);
    k_31 = F_tr1((t(i)+0.5*h),(r(1,i)+0.5*h*k_21),(r(2,i)+0.5*h*k_22));
    k_32 = F_tr2((t(i)+0.5*h),(r(1,i)+0.5*h*k_21),(r(2,i)+0.5*h*k_22));
    k_41 = F_tr1((t(i)+h),(r(1,i)+k_31*h),(r(2,i)+k_32*h));
    k_42 = F_tr2((t(i)+h),(r(1,i)+k_31*h),(r(2,i)+k_32*h));
    %First ODE result
    r(1,i+1) = r(1,i) + (1/6)*(k_11+2*k_21+2*k_31+k_41)*h;  % main equation RK4
    %Second ODE result
    r(2,i+1) = r(2,i) + (1/6)*(k_12+2*k_22+2*k_32+k_42)*h;  % main equation RK4
end

end
