clc;
clear all;
close all;
NoB = 1;
v = [0.01 0.1 0.3];
K = -5:1:5;
for V = 1:length(v)
for k = 1:length(K)
syms y a hs hp a1;
P = a/hs; PP = a1/hs;
fhp = (1/(1-v(V)))*exp(-hp/(1-v(V)));
fhs = exp(-hs);
Pin1 = int(fhs,hs,(hp+v(V))/y,inf);
Pin = int(fhp*Pin1,hp,0,inf) == 0.8;
S = solve(Pin, y); Y(V) = double(S);
P1 = P*(hp+v(V))*fhp;
Q0 = int(P1,hp,0,Y(V)*hs - v(V));
Q = int(fhs*Q0, hs , v(V)/Y(V),inf) == 10^(K(k)/10);
S1 = solve(Q, a); A(V) = double(S1);
Cout(k,V) = log(1+(A(V)/NoB))*((Y(V)*exp(-v(V)/Y(V)))/(Y(V)+1-v(V)));

[Y1 Ctifr(k,V)] = fminbnd(@(Y1) -(log(1+((10^(K(k)/10))/(-(Y1*(1-v(V))/(Y1+1-v(V)))*exp(-v(V)/Y1) - ei(-v(V)/Y1) + (1-v(V))*exp(v(V)/(1-v(V)))*ei((-v(V)/(1-v(V))) - v(V)/Y1))))*((Y1*exp(-v(V)/Y1))/(Y1+1-v(V)))),0,10000);

end
end
plot(K,(Cout(:,1)),'-o');
hold on;
plot(K,-Ctifr(:,1),'-*');
hold on;
plot(K,(Cout(:,2)),'-o');
hold on
plot(K,-Ctifr(:,2),'-*');
hold on
plot(K,(Cout(:,3)),'-o');
hold on
plot(K,-Ctifr(:,3),'-*');
legend({'var = 0.01 cout','var = 0.01 ctifr','var = 0.1 cout','var = 0.1 ctifr','var = 0.3 cout','var = 0.3 ctifr'},'Location','northwest');
xlabel('Qavg (dB)');
ylabel('Ctifr/B');
title('Ctifr/B vs Qavg - average received power constraint');