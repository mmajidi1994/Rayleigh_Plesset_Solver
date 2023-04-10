function [T,Tz,Tzz] = initial(n,z,t_inf)

syms TT [1,n-1] 

syms k x 
TT=chebyshevT(2*(1:n-1)-1, x);
temp=zeros(n-1,n-1);
% Initial condition
for i=1:n-1
   for j=1:n-1
    temp(i,j)=chebyshevT(2*j-1, z(i+1)); %...
  %-subs(diff(chebyshevT(2*j-1, x),x), 1)/subs(diff(chebyshevT(2*(n-1)-1, x),x), 1);
   end
end

T0=t_inf*ones(n-1,1);
a=inv(temp)*T0;

T(x)= dot(a,TT);
Tz(x)=diff(T,x);
Tzz(x)=diff(T,x,x);
%plot(z(2:n),T(z(2:n)))

end
