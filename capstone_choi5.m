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

syms Ts z
t_cond=2400; HAs(1)=0; Tss(1)=250;
for i=2:366
    del_z=h(i-1)/10;
   
    equ2=ki*(Ta(i+30)-Ts)/del_z;
   
    qsat=0.6+0.06*(Ts-273);
    SHF=rhoa*cp*CD*U(i+30)*(Ts-Ta(i+30));
    equ3=SHF+eps*sigm*Ts^4-(1-alp)*(1-I0)*fswdn(i+30)-flwdn(i+30);
   
    equ=equ3-equ2;
    Ts_sl=vpasolve(equ,Ts,[220 300]);
    Tss(i)=double(Ts_sl);
    HAs(i)=subs(equ3,Tss(i));
   
    if HAs(i)<=0 | Tss(i)>=273
        Tss(i)=Tss(i-1); HAs(i)=0;
    end
   
    equ1io=I0*(1-alp)*fswdn(i+30)*exp(-kap*z);
    equ1i=diff(equ1io,z);
    RH=subs(equ1i,z,del_z);
   
    Tm(1,i)=Tss(i-1);
    for j=2:10
        Tm(j,i)=Tm(j,i-1)+(t_cond/(rhoi*cpi))*...
            (ki*(Tm(j+1,i-1)-2*Tm(j,i-1)+Tm(j-1,i-1))/del_z^2-RH);
    end
   
    T=Tm(:,i)';
    h(i)=h(i-1)+del_t/Lfib*(ki*(T(11)-T(10))/del_z-Fw-I0*(1-alp)...
        *fswdn(i+30)*exp(kap*h(i-1)));
end
