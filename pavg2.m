clc;
clear all;
close all;
NoB = 1;
v = [0.01 0.05 0.1];
K = -5:1:5;

for j = 1:length(v)
for i = 1:length(K)
    syms a hs hp y;
    P = a./hs;
    fhp = (1/(1-v(j)))*exp(-hp/(1-v(j)));
    fhs = exp(-hs);
    
    Pin1 = int(fhs,hs,(hp+v(j))/y,inf);
    Pin = int(Pin1*fhp,hp,0,inf) == 0.8;
    Ya = solve(Pin,y); Y(j) = double(Ya);
    
    Q0 = int(P*(hp+v(j))*fhp,hp,0,Y(j)*hs - v(j));
    Q = int(Q0*fhs,hs,v(j)/Y(j),inf) == 10^(K(i)/10);
    aa = solve(Q,a); A(j,i) = double(aa);
    
    PP = (A(j,i)/hs)*fhp*fhs;
    Pavg1 = int(PP,hp,0,Y(j)*hs-v(j));
    Pavg2 = int(Pavg1,hs,v(j)/Y(j),inf);
    Pavg(j,i) = double(Pavg2);
    
end
end

plot(K,10*log10(Pavg(1,:)));
hold on;
plot(K,10*log10(Pavg(2,:)))
hold on;
plot(K,10*log10(Pavg(3,:)))
axis([-5,5,-10,15])