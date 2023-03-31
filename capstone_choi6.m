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
rhoi = 913;            
rhoa = 1.275;          
sigma = 5.67*10^-8;   
k_i = 2;              
I0 = 0.17;                          
alpha = 0.5;                  
kappa = 1.5;           
Fw = 2; 
k_s = 0.31;

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

day_hi = 2.*ones(10,365);
delta_t = 3600;
day_T = ones(10,365,11); z = ones(10,365,11); delta_z = ones(10,365); day_HA = ones(10,365); SHF_c = ones(10,365);
T_n = T; z0_n = 0:0.2:2; dz_n = 0.2; hi_n = 2;

for iyear = 1:10
    for ihour = 1:365*24
        z0 = z0_n; dz = dz_n; T = T_n; hi = hi_n;
        
        syms Ts

        eq2 = k_i*((T(2)-Ts)/dz);

        qsat = 0.6 + 0.06*(Ts-273);
        SHF = rhoa*cp*Cd*U(ihour)*(Ts-Ta(ihour));
        LHF = 0;

        eq3 = sigma*Ts^4 - flwdn(ihour) - (1-alpha)*(1-I0)*fswdn(ihour) + SHF + LHF;

        sol_Ts = vpasolve(eq3==eq2, Ts, [220 300]);
        T(1) = double(sol_Ts);
        if T(1) >= 273
            T(1) = 273;
        end
        HA = double(subs(eq3, Ts, T(1)));
        
        T_n(1) = T(1); T_n(11) = 271;
        for i = 2:10
            T_n(i) = T(i) + delta_t/(c_pi*rhoi)*(k_i*((T(i+1)-2*T(i)+T(i-1))/dz^2) ...
                - I0*(1-alpha)*fswdn(ihour)*(-kappa)*exp(-kappa*z0(i)));            
            if T_n(i) >= 273
                T_n(i) = 273;
            end
        end
        
        fs = k_s*(T(1)-Ts)/hs);
        
        if HA < 0
            hi_n = hi + HA*delta_t/Lfi;
        else
            hi_n = hi + delta_t/Lfib*(k_i*(T_n(end)-T_n(end-1))/dz ...
                - Fw - I0*(1-alpha)*fswdn(ihour)*exp(-kappa*hi));
        end
        dz_n = hi_n/10;
        z0_n = 0:dz_n:hi_n;
        
        if rem(ihour,24) == 1
            day_HA(iyear,round(ihour/24)+1) = HA;
            SHF_c(iyear,round(ihour/24)+1) = double(subs(SHF, Ts, T_n(1)));
            day_T(iyear,round(ihour/24)+1,:) = T_n;
            z(iyear,round(ihour/24)+1,:) = z0_n;
            delta_z(iyear,round(ihour/24)+1) = dz_n;
            day_hi(iyear,round(ihour/24)+1) = hi_n;
        end
    end
end

%% plot
figure(1)
set(gcf, 'position', [300 300 800 400])
subplot(1,2,1)
plot(1:365, day_T(10,:,1))
xlabel('day')
title('Ts (K)')

subplot(1,2,2)
plot(1:365,day_hi(10,:))
xlabel('day')
title('ice thickness (m)')

figure(2)
set(gcf, 'position', [500 300 800 400])
subplot(1,3,1)
plot(squeeze(day_T(10,(31+14),:)),squeeze(z(10,(31+14),:)),'-o','linewidth',3)
xlabel('Temp (K)')
ylabel('depth (m)')
title('2월 14일')
ylim([0 2.3])
axis ij

subplot(1,3,2)
plot(squeeze(day_T(10,(232+5),:)),squeeze(z(10,(232+5),:)),'-o','linewidth',3)
xlabel('Temp (K)')
ylabel('depth (m)')
title('8월 5일')
ylim([0 2.3])
axis ij

subplot(1,3,3)
plot(squeeze(day_T(10,(273+25),:)),squeeze(z(10,(273+25),:)),'-o','linewidth',3)
xlabel('Temp (K)')
ylabel('depth (m)')
title('10월 25일')
ylim([0 2.3])
axis ij















