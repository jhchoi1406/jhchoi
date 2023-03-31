clear; clc;

load('force1D_origin.mat')

% method 1 : added extrapolation %

month = 1:12;
day0 = 0.5394:0.0329:12.5264; % 2월은 28일로 가정 , 남은 일수 30일

value01 = interp1(month,avg_flwdn,day0,'spline','extrap');
value02 = interp1(month,avg_fswdn,day0,'spline','extrap');
value03 = interp1(month,avg_Qref,day0,'spline','extrap');
value04 = interp1(month,avg_Ta,day0,'spline','extrap');
value05 = interp1(month,avg_U,day0,'spline','extrap');

figure(1)
subplot(3,2,1)
plot(month,avg_flwdn,'r>',day0,value01,'b--')
subplot(3,2,2)
plot(month,avg_fswdn,'r>',day0,value02,'b--')
subplot(3,2,3)
plot(month,avg_Qref,'r>',day0,value03,'b--')
subplot(3,2,4)
plot(month,avg_Ta,'r>',day0,value04,'b--')
subplot(3,2,5)
plot(month,avg_U,'r>',day0,value05,'b--')

title('method1')

% method 2 : used interpolation only %

day = 1:0.0329:12;

value1 = interp1(month,avg_flwdn,day);
value2 = interp1(month,avg_fswdn,day);
value3 = interp1(month,avg_Qref,day);
value4 = interp1(month,avg_Ta,day);
value5 = interp1(month,avg_U,day);

m = [1, 2];
d = linspace(1,2,32);

v1 = [avg_flwdn(end), avg_flwdn(1)];
p1 = interp1(m,v1,d);
t1 = [p1(18:end-1) value1 p1(2:17)];
v2 = [avg_fswdn(end), avg_fswdn(1)];
p2 = interp1(m,v2,d);
t2 = [p2(18:end-1) value2 p2(2:17)];
v3 = [avg_Qref(end), avg_Qref(1)];
p3 = interp1(m,v3,d);
t3 = [p3(18:end-1) value3 p3(2:17)];
v4 = [avg_Ta(end), avg_Ta(1)];
p4 = interp1(m,v4,d);
t4 = [p4(18:end-1) value4 p4(2:17)];
v5 = [avg_U(end), avg_U(1)];
p5 = interp1(m,v5,d);
t5 = [p5(18:end-1) value5 p5(2:17)];

figure(2)
subplot(3,2,1)
plot(month,avg_flwdn,'r>',day0,t1,'b--')
subplot(3,2,2)
plot(month,avg_fswdn,'r>',day0,t2,'b--')
subplot(3,2,3)
plot(month,avg_Qref,'r>',day0,t3,'b--')
subplot(3,2,4)
plot(month,avg_Ta,'r>',day0,t4,'b--')
subplot(3,2,5)
plot(month,avg_U,'r>',day0,t5,'b--')

figure(3)
plot(month,avg_flwdn,'r>',day0,value01,'b')
grid on; hold on;
plot(month,avg_flwdn,day0,t1,'k')
hold off;

legend('mean','method1','method2')








