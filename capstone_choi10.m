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
snow=0.01*interp1(snt,sn0,dayt,'spline');

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
k_i = 2;              
i0 = 0.17;                          
alpha = 0.5;                  
kappa = 1.5;           
Fw = 2; 
k_s = 0.31;
z0 = 0:0.2:2;
hi0 = 2;
delta_t = 3600;
hi = 1.38;
t_cond = 2400;

T = [250 260 265 270 275 280 288 283 279 277 271];

hs(1)=0; hi(1)=hi; T(1)=0;

for i = 1:365*24

    if hs==0
        syms Ts          
        qsat = 0.6 + 0.06*(Ts-273);
        SHF = rhoa*cp*Cd*U(i)*(Ts-Ta(i));
        LHF = 0;
            
        eq2 = ki*((T(2)-Ts)/dz);
        eq3 = sigma*Ts^4 - flwdn(ihour) - (1-alpha_i)*(1-I0)*fswdn(ihour) + SHF + LHF;

            % 표층온도 T(1) 구하기
            sol_Ts = vpasolve(eq3-eq2, Ts, [220 300]);
            Tsnow = nan;
            T(1) = double(sol_Ts);
            if T(1) >= 273
                T(1) = 273;
            end
            HA = double(subs(eq3, Ts, T(1)));
            
            T_n(1) = T(1);
    syms Ts z
    
    delta_z=hi(i-1)/10;
    hsns=hs(i-1)+snow(i-1);
    eq1=k_i*k_s*(Tf-Ts)/(k_i*hsns+k_s*hi(i-1));
    
    SHF = rhoa*cp*Cd*U(i-1)*(Ts-Ta(i-1));
    LHF = 0;

    eq2 = sigma*Ts^4 -fswdn(i-1) - flwdn(i-1) + SHF + LHF ;
    eq3 = eq1-eq2;
    Ts0 = vpasolve(eq3, Ts, [220 300]);
    Ts = double(Ts0);
    HA = subs(eq1, Ts);

    hs(i)=hs(i-1)+HA*delta_t/Lfs;
    if hs(i)<=0.001
        hs(i)=0;
        if HA<0 | Ts>=273
           T1=T(1); HA=0;
        end
    else
        T1=HA*hs(i)/k_s+Ts;
    end
    T(1)=T1;

    equ1io=i0*(1-alpha)*fswdn(i-1)*exp(-kappa*z);
    equ1i=diff(equ1io,z);
    RH=subs(equ1i,z,delta_z);

   
    for j = 2:10
        T(j) = T(j-1) + t_cond/(c_pi*rhoi)*(k_i*((T(j+1)-2*T(j)+T(j-1))/delta_z^2-RH)); 
    end
    
    hi(i)=hi(i-1)+delta_t/Lfib*(k_i*(T(11)-T(10))/delta_z-Fw-i0*(1-alpha)...
        *fswdn(i-1)*exp(-kappa*hi(i-1)));

end

plot(dayt,hi(2:8761),dayt,hi(2:8761)+hs(2:8761))











