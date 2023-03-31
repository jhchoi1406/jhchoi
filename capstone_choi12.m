clear;

load('force1D_origin.mat')
sn = open('snowfall_1.mat')

avg_flwdn = [avg_flwdn; avg_flwdn];
avg_fswdn = [avg_fswdn; avg_fswdn];
avg_Qref = [avg_Qref; avg_Qref];
avg_Ta = [avg_Ta; avg_Ta];
avg_U = [avg_U; avg_U];

xq = 0:24/(365*24*2):24;
flwdn0 = interp1(avg_flwdn, xq);
fswdn0 = interp1(avg_fswdn, xq);
Qref0 = interp1(avg_Qref, xq);
Ta0 = interp1(avg_Ta, xq);
U0 = interp1(avg_U, xq);

idx = find(~isnan(flwdn0),1);
flwdn0(1:idx-1) = flwdn0(365*24+1:365*24+idx-1);
fswdn0(1:idx-1) = fswdn0(365*24+1:365*24+idx-1);
Qref0(1:idx-1) = Qref0(365*24+1:365*24+idx-1);
Ta0(1:idx-1) = Ta0(365*24+1:365*24+idx-1);
U0(1:idx-1) = U0(365*24+1:365*24+idx-1);

flwdn = flwdn0(1:365*24);
fswdn = fswdn0(1:365*24);
Qref = Qref0(1:365*24);
Ta = Ta0(1:365*24);
U = U0(1:365*24);

sn0=[sn.snow_interp]/24';
snt=24*(1:365)-12;
dayt=1:8760;
snow=interp1(snt,sn0,dayt,'spline');

cp = 1005;              
c_pi = 2092;           
Cd = 0.0013;                  
Ls = 2.265*10^6;       
Lfi = 3.014*10^8;      
Lfib = 2.679*10^8;  
Lfs = 1.097*10^8;
rhoi = 913;            
rhoa = 1.275;    
Tf = 271.3;
sigma = 5.67*10^-8;   
ki = 2;              
i0 = 0.17;                          
alpha = 0.5;  
alpha_s = 0.8;
kappa = 1.5;           
Fw = 2; 
ks = 0.31;
z0 = 0:0.2:2;
hi0 = 2;
delta_t = 3600;
hi = 1.38;
t_cond = 2400;

T = [250 260 265 270 275 280 288 283 279 277 271];

day_T = ones(10,365,11); day_Tsnow = ones(10,365); day_hi = ones(10,365); day_hs = zeros(10,365);
z = ones(10,365,11); delta_z = ones(10,365); day_HA = ones(10,365); SHF_c = ones(10,365);

z0=linspace(0,1,11); z=ones(1,365,11); T0=T; hi0=1; dz0=hi0/3; hs0=0;

for i = 1:365*24 
    zt = z0; dz = dz0; T = T0; hi = hi0; hs = hs0;

    if hs == 0 
        syms Ts          
        qsat = 0.6 + 0.06*(Ts-273);
        SHF = rhoa*cp*Cd*U(i)*(Ts-Ta(i));
        LHF = 0;
       
        eq1 = ki*((T(2)-Ts)/dz);
        eq2 = sigma*Ts^4 - flwdn(i) - (1-alpha)*(1-i0)*fswdn(i) + SHF + LHF;

        sol_Ts = vpasolve(eq2-eq1, Ts, [220 300]);
        Tsnow = nan;

        T(1) = double(sol_Ts);
        if T(1) >= 273
            T(1) = 273;
        end
        HA = double(subs(eq2, Ts, T(1)));
            
        T0(1) = T(1);
            
        T0(11) = 271;    
        for j = 2:10
            T0(j) = T(j) + delta_t/(c_pi*rhoi)*(ki*((T(j+1)-2*T(j)+T(j-1))/dz^2) ...
                - i0*(1-alpha)*fswdn(i)*(-kappa)*exp(-kappa*zt(j)));            
            if T0(j) >= 273
                T0(j) = 273;
            end
        end                    

        HAs = HA;
        if HAs > 0
            HAs = 0;
            hs0 = hs + snow(i) + HAs*delta_t/Lfs;
            hi0 = hi + delta_t/Lfib*(ki*(T0(end)-T0(end-1))/dz ...
                    -Fw-i0*(1-alpha)*fswdn(i)*exp(-kappa*hi));
        else
            hs0 = hs + snow(i) + HAs*delta_t/Lfs;
            if hs0 <= 0
                HA1 = -(hs+snow(i))*Lfs/delta_t;
                HA2 = HAs - HA1;
                hs0= 0;
                hi0 = hi + HA2*delta_t/Lfi;
            end    
        end
        dz0 = hi0/10;
        z0 = 0:dz0:hi0;
            
    else
        syms Ts          
        qsat = 0.6 + 0.06*(Ts-273);
        SHF = rhoa*cp*Cd*U(i)*(Ts-Ta(i));
        LHF = 0;
            
        eq1 = ks*((T(1)-Ts)/hs);
        eq2 = sigma*Ts^4 - flwdn(i) - (1-alpha_s)*(1-i0)*fswdn(i) + SHF + LHF;

        sol_Ts = vpasolve(eq2-eq1, Ts, [220 300]);
        Tsnow = double(sol_Ts);
        if Tsnow >= 273
            Tsnow = 273;
        end
        HA = double(subs(eq2, Ts, Tsnow));

        T0(1) = (hs*ki*T(2)+dz*ks*Tsnow)/(hs*ki+dz*ks);
            
        T0(11) = 271;    
        for j = 2:10
            T0(j) = T(j) + delta_t/(c_pi*rhoi)*(ki*((T(j+1)-2*T(j)+T(j-1))/dz^2) ...
                - i0*(1-alpha)*fswdn(i)*(-kappa)*exp(-kappa*z0(j)));            
            if T0(j) >= 273
                T0(j) = 273;
            end
        end                    

        HAs = HA;
        if HAs > 0 
            HAs = 0;
            hs0 = hs + snow(i) + HAs*delta_t/Lfs;
            hi0 = hi + delta_t/Lfib*(ki*(T0(end)-T0(end-1))/dz ...
                    - Fw - i0*(1-alpha_s)*fswdn(i)*exp(-kappa*hi));
        else 
            hs0 = hs + snow(i) + HAs*delta_t/Lfs;
            if hs0 <= 0
                HA1 = -(hs+snow(i))*Lfs/delta_t;
                HA2 = HAs - HA1;
                hs0= 0;
                hi0 = hi + HA2*delta_t/Lfi;
            end    
        end
        dz0 = hi0/10;
        z0 = 0:dz0:hi0;
    end

    end











