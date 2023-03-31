clear;

load('historical_forcing+snow(daily).mat')

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
delta_t = 86400;

T = [250 260 271];

z0_1 = linspace(0,1,3);
% hi0=0.1; dz0=hi0/2; hs0=0;
hi=0.1; dz=hi/2; hs=0;
HA_t = []; hi_t = []; hs_t = [];
Ts01 = [];
for i = 1:365
%     z0 = z0_1; dz = dz0; T0 = T; hi = hi0; hs = hs0;
    
    if hs==0 % 눈x
        
        
        T(3) = 271;
        for k = 2
            T(k) = T(k) + delta_t/(c_pi*rhoi)*(ki*((T(k+1)-2*T(k)+T(k-1))/dz^2) ...
                - i0*(1-alpha)*fswdn(i)*(-kappa)*exp(-kappa*z0(k)));            
            if T(k) >= 273
                T(k) = 273;
            end
        end
        
        
        syms Ts
        qsat = 0.6 + 0.06*(Ts-273);
        SHF = rhoa*cp*Cd*U(i)*(Ts-Ta(i));
        LHF = 0;
        
        eq1 = ki*((T(2)-Ts)/dz);
        eq2 = sigma*Ts^4 - flwdn(i) - (1-alpha)*(1-i0)*fswdn(i) + SHF + LHF;

        Ts0 =abs(vpasolve(eq2==eq1, Ts));
        if Ts0 >=273
            Ts0 = 273
        else
            Ts0 = Ts0(find(Ts0>=220 & Ts0<=300));
        end

        Ts01 = [Ts01; Ts0];
%         HA = double(subs(eq2, Ts, T(1)));
        if HA > 0
                HA = 0;
                hs0 = hs + snow(i) + HA*delta_t/Lfs;
                hi0 = hi + delta_t/Lfib*(ki*(T0(end)-T0(end-1))/dz ...
                        - Fw - i0*(1-alpha)*fswdn(i)*exp(-kappa*hi));
            else 
                hs0 = hs + snow(i) + HA*delta_t/Lfs;
                if hs0 <= 0
                    HA1 = -(hs+snow(i))*Lfs/delta_t;
                    HA2 = HA - HA1;
                    hs0= 0;
                    hi0 = hi + HA2*delta_t/Lfi;
                end    
            end
            dz0 = hi0/2;
            z0 = 0:dz0:hi0;
        else
            syms Ts
            qsat = 0.6 + 0.06*(Ts-273);
            SHF = rhoa*cp*Cd*U(i)*(Ts-Ta(i));
            LHF = 0;

            eq1 = ks*((T(1)-Ts)/hs);
            eq2 = sigma*Ts^4 - flwdn(i) - (1-alpha_s)*(1-i0)*fswdn(i) + SHF + LHF;

            Ts0 = vpasolve(eq2-eq1, Ts, [220 300]);
            Tws = double(Ts0);
            if Tws >= 273
                Tws = 273;
            end
            HA = double(subs(eq2, Ts, Tws));

            T0(1) = (hs*ki*T(2)+dz*ks*Tws)/(hs*ki+dz*ks);

            T0(3) = 271;  
            for k = 2
                T0(k) = T(k) + delta_t/(c_pi*rhoi)*(ki*((T(k+1)-2*T(k)+T(k-1))/dz^2) ...
                    - i0*(1-alpha_s)*fswdn(i)*(-kappa)*exp(-kappa*z0(k)));            
                if T0(k) >= 273
                    T0(k) = 273;
                end
            end

            if HA > 0
                HA = 0;
                hs0 = hs + snow(i) + HA*delta_t/Lfs;
                hi0 = hi + delta_t/Lfib*(ki*(T0(end)-T0(end-1))/dz ...
                        - Fw - i0*(1-alpha_s)*fswdn(i)*exp(-kappa*hi));
            else 
                hs0 = hs + snow(i) + HA*delta_t/Lfs;
                if hs0 <= 0
                    HA1 = -(hs+snow(i))*Lfs/delta_t;
                    HA2 = HA - HA1;
                    hs0= 0;
                    hi0 = hi + HA2*delta_t/Lfi;
                end    
            end
            dz0 = hi0/2;
            z0 = 0:dz0:hi0;

        end
        HA_t = [HA_t; HA];
        hi_t = [hi_t; hi0];
        hs_t = [hs_t; hs0];

        
        
    else
        
    end
end















































