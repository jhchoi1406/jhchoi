clear ;

z = 0:0.2:2 ;
T = [250 251 252 254 256 259 262 265 267 269 271] ;
rhoi = 913 ;         
cpi = 2093 ;
ki = 2.22 ;
delz = 0.2 ;
delt = 1000 ;

n = 1;
while n ~= 10000000
    for i = 2:length(T)-1
        T(i) = T(i) + delt*((ki/(cpi*rhoi))*((T(i+1)-2*T(i)+T(i-1))/delz^2));
        T(1) = 250;
        T(end) = 271;
        n = n+1;
    end
end
