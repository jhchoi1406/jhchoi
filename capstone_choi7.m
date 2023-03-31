clear;

load('force1D_origin.mat')
load('snowfall_1.mat')

month = 1:12;
day = 1:365;
day0 = 0:12/(24*365):12;

flwdn0 = interp1(month,avg_flwdn,day0,'spline','extrap');
flwdn = flwdn0(2:end);
fswdn0 = interp1(month,avg_fswdn,day0,'spline','extrap');
fswdn = fswdn0(2:end);
Qref0 = interp1(month,avg_Qref,day0,'spline','extrap');
Qref = Qref0(2:end);
Ta0 = interp1(month,avg_Ta,day0,'spline','extrap');
Ta = Ta0(2:end);
U0 = interp1(month,avg_U,day0,'spline','extrap');
U = U0(2:end);
snow0 = interp1(day,snow_interp,day0, 'spline','extrap');
snow = snow0(2:end);

cp = 1005;              
c_pi = 2092;           
Cd = 0.0013;                  
Ls = 2.265*10^6;       
Lfi = 3.014*10^8;      
Lfib = 2.679*10^8;  
Lfs = 1.097*10^8;
rhoi = 913;            
rhoa = 1.275;          
sigma = 5.67*10^-8;   
k_i = 2;              
i0 = 0.17;                          
alpha = 0.5;                  
kappa = 1.5;           
Fw = 2; 
k_s = 0.31;
z0 = 0:0.2:2;
hi0 = 2;
delta_t = 3600;
hs0 = 1.38;

dz = 0.2;
dt = 1000;
T = [250 260 265 270 275 280 288 283 279 277 271];
n = 1; 
while n ~= 5000
    for i = 2:10
        T(i) = T(i) + dt*((k_i/(c_pi*rhoi))*((T(i+1)-2*T(i)+T(i-1))/dz^2));
        T(1) = 250;
        T(end) = 271;
    end
    n = n+1;
end

hi_t = [];
for i = 1:365*24
    syms Ts
    
    eq1 = k_s*(T(1)-Ts)/dz ;
    
    SHF = rhoa*cp*Cd*U(i)*(Ts-Ta(i));
    LHF = 0;

    eq2 = sigma*Ts^4 - (1-alpha)*(1-i0)*fswdn(i) - flwdn(i) + SHF + LHF ;
    Ts0 = vpasolve(eq1 == eq2, Ts, [220 300]);
    T(1) = double(Ts0);
    if T(1) >= 273
       T(1) = 273;
    end
    HA = double(subs(eq2, Ts, T(1)));
    
    for j = 2:10
        T(j) = T(j) + delta_t/(c_pi*rhoi)*(k_i*((T(j+1)-2*T(j)+T(j-1))/dz^2) ...
            - i0*(1-alpha)*fswdn(i)*(-kappa)*exp(-kappa*z0(j)));    
        if T(j) >=273
               T(j) = 273;
        end
    end
    T(11) = 271;
    
    if HA < 0
        hi = hi0 + HA*delta_t/Lfi;
    else
        hi = hi0 + delta_t/Lfib*(k_i*(T(end)-T(end-1))/dz ...
         - Fw - i0*(1-alpha)*fswdn(i)*exp(-kappa*hi0));
    end
    hi_t = [hi_t; hi];

end

figure(1)
plot(hi_t)
xlabel('hour')
ylabel('ice thickness')





















