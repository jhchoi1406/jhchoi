 +clear;

load('force1D_origin.mat')
load('snowfall_1.mat')

avg_flwdn = [avg_flwdn; avg_flwdn];
avg_fswdn = [avg_fswdn; avg_fswdn];
avg_Qref = [avg_Qref; avg_Qref];
avg_Ta = [avg_Ta; avg_Ta];
avg_U = [avg_U; avg_U];

xq = 0:24/(365*2):24;
flwdn0 = interp1(avg_flwdn, xq);
fswdn0 = interp1(avg_fswdn, xq);
Qref0 = interp1(avg_Qref, xq);
Ta0 = interp1(avg_Ta, xq);
U0 = interp1(avg_U, xq);

idx = find(~isnan(flwdn0),1);
flwdn0(1:idx-1) = flwdn0(365+1:365+idx-1);
fswdn0(1:idx-1) = fswdn0(365+1:365+idx-1);
Qref0(1:idx-1) = Qref0(365+1:365+idx-1);
Ta0(1:idx-1) = Ta0(365+1:365+idx-1);
U0(1:idx-1) = U0(365+1:365+idx-1);

flwdn = flwdn0(1:365);
fswdn = fswdn0(1:365);
Qref = Qref0(1:365);
Ta = Ta0(1:365);
U = U0(1:365);

snow = 0.01*snow_interp;

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
        

        
        
    else
        
    end
end














































