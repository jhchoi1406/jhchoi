
clear;

cp = 1005;              
c_pi = 2092;           
Cd = 0.0013;                            
Ls = 2.265*10^6;       
Lfi = 3.014*10^8;      
rhoi = 913;            
rhoa = 1.275;          
sigma = 5.67*10^-8;    
U = 2.56;              
Ta = 249.85;            
qa = 0.57;              
F_LW = 182.1;          
k_i = 2;               

T = [250 260 265 270 275 280 288 283 279 277 271];
z = 0:0.2:2; 
delta_z = 0.2;
delta_t = 1000;

n = 1; nn = 1; t = 0; hi = 1;
while n ~= 10
    syms Ts HA

    eq2 = k_i*((T(2)-Ts)/delta_z) == HA;

    qsat = 0.6 + 0.06*(Ts-273);
    SHF = rhoa*cp*Cd*U*(Ts-Ta);
    LHF = rhoa*Ls*Cd*U*(qsat*Ts-qa);
    eq3 = sigma*Ts^4 - F_LW + SHF + LHF == HA;

    [sol_Ts, sol_HA] = vpasolve([eq2 eq3],[Ts HA]);
    T(1) = sol_Ts(3);

    while nn ~= 5000
        for i = 2:length(T)-1
            T(i) = T(i) + delta_t*((k_i/(c_pi*rhoi))*((T(i+1)-2*T(i)+T(i-1))/delta_z^2));
            T(end) = 271;
        end
        nn = nn+1;
    end
    n = n+1;
end

for j = 1:30
    hi(j+1) = sol_HA(3)*delta_t/Lfi + hi(j);
    t(j+1) = t(j) + delta_t;
end

figure(1)
plot(T,z)
axis ij
xlabel('Temperature')
ylabel('depth')
title('T')

figure(2)
plot(t,hi)
xlabel('time')
ylabel('ice thickness')
title('hi')
