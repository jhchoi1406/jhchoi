clear;

load('historical_forcing+snow(daily).mat')

fswdn = fswdn_interp; flwdn = flwdn_interp; U = U_interp; Ta = Ta_interp; Qref = Qref_interp; snow = snow_interp;

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
delta_t = 86400;

T = [250 260 271];
hi_n=0.1; dz_n=hi_n/2; hs_n=0;
hi_t = []; hs_t = []; HA_t = []; T_t = []; T_n = T; dz_t=[];

for i = 1:365
    dz = dz_n; hi = hi_n; hs = hs_n; T = T_n;
    if hs == 0
        T_n(3) = 271;
        for k = 2
            if hi>=0.4
                T_n(k) = T(k) + delta_t/(c_pi*rhoi)*(ki*((T(k+1)-2*T(k)+T(k-1))/dz^2) ...
                    - i0*(1-alpha)*fswdn(i)*(-kappa)*exp(-kappa*hi));      
                if T_n(k) >= 273
                    T_n(k) = 273;
                end
            elseif hi<0.4 & hi>=0.01
                T_n(k) = 
            else
                T_n(:) = 273;
            end
        end

        syms Ts
        qsat = 0.6 + 0.06*(Ts-273);
        SHF = rhoa*cp*Cd*U(i)*(Ts-Ta(i));
        LHF = 0;
        
        eq1 = ki*((T_n(2)-Ts)/dz);
        eq2 = sigma*Ts^4 - flwdn(i) - (1-alpha)*(1-i0)*fswdn(i) + SHF + LHF;

        Ts=double(vpasolve(eq2-eq1, Ts, [220 300]));
        T_n(1) = Ts;
        Tsnow = T_n(1);
        syms Ts
        HA = double(subs(eq2, Ts, T_n(1)));

        if HA>0
            HA = 0;
            hs_n = hs + snow(i) + HA*delta_t/Lfs;
            hi_n = hi + delta_t/Lfib*(ki*(T_n(end)-T_n(end-1))/dz-Fw-i0*(1-alpha_s)*fswdn(i)*exp(-kappa*hi));

        else
            hs_n = hs + snow(i) + HA*delta_t/Lfs;
            if hs_n <=0
                hs_n = 0;
                hi_n = hi + (HA+(hs+snow(i))*Lfs/delta_t)*delta_t/Lfi;
            else 
                hi_n = hi;
            end
        end
        dz_n = hi_n/2;

    else
        T_n(3) = 271;
        for k = 2
            T_n(k) = T(k) + delta_t/(c_pi*rhoi)*(ki*((T(k+1)-2*T(k)+T(k-1))/dz^2) ...
                - i0*(1-alpha)*fswdn(i)*(-kappa)*exp(-kappa*hi));           
            if T_n(k) >= 273
                T_n(k) = 273;
            end
        end

        syms Ts
        qsat = 0.6 + 0.06*(Ts-273);
        SHF = rhoa*cp*Cd*U(i)*(Ts-Ta(i));
        LHF = 0;

        eq1 = ki*ks*(Tf-Ts)/(ki*hs+ks*hi);
        eq2 = sigma*Ts^4 - flwdn(i) - (1-alpha_s)*(1-i0)*fswdn(i) + SHF + LHF;

        Ts=double(vpasolve(eq2-eq1, Ts, [220 300]));

        Tsnow = Ts;
        if Tsnow >= 273
            Tsnow = 273;
        end
        syms Ts
        HA = double(subs(eq2, Ts, Tsnow));

        T_n(1) = (hs*ki*T_n(2)+dz*ks*Tsnow)/(hs*ki+dz*ks);
        i

        if HA>0
            HA = 0;
            hs_n = hs + snow(i) + HA*delta_t/Lfs;
            hi_n = hi + delta_t/Lfib*(ki*(T_n(end)-T_n(end-1))/dz-Fw-i0*(1-alpha_s)*fswdn(i)*exp(-kappa*hi));

        else
            hs_n = hs + snow(i) + HA*delta_t/Lfs;
            if hs_n <=0
                hs_n = 0;
                hi_n = hi + (HA+(hs+snow(i))*Lfs/delta_t)*delta_t/Lfi;
            else 
                hi_n = hi;
            end
        end
        dz_n = hi_n/2;
    end
    hi_t = [hi_t; hi_n];
    hs_t = [hs_t; hs_n];
    HA_t = [HA_t; HA];
    T_t = [T_t; T_n];
    dz_t = [dz_t; dz_n];
   
end
figure(1)
plot(hi_t)
figure(2)
plot(hs_t)
