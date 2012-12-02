%load out/rho0000.dat;
newplot;
hold on;
load out/rho0001.dat;
plot(rho0001(100,:),'o-');
load out1/rho0001.dat;
plot(rho0001(50,:),'o-r');
hold off;
