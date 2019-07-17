clc;
clear all;
close all;
NoB = 1;
v = [0,0.01,0.05,0.1];
K = -5:1:5;
for V = 1:length(v)
for k = 1:length(K)
syms y hp hs;
P = (y/(hp+v(V))) - (NoB/hs);
fhp = (1/(1-v(V)))*exp(-hp/(1-v(V)));
fhs = exp(-hs);
P1 = P*(hp+v(V))*fhp*fhs;
Q0 = int(P1,hp,0,(y*hs/NoB)-v(V));
Q = int(Q0,hs,NoB*v(V)/y,inf) == 10^(K(k)/10);
S = solve(Q,y); Si(k) = double(S);

for i = 1:1000
    hpc = exprnd(1-v(V));
    if hpc >= 0
        hpc = hpc;
    else hpc = 0;
    end
    hsc = exprnd(1);
    
    p(i) = (Si(k)/(hpc+v(V))) - (NoB/hsc);
    if p(i) >= 0
        p(i) = p(i);
    else p(i) = 0; end
    
end
pav(k,V) = mean(p);

pavt(k,V) = -((Si(k)/(1-v(V)))+NoB)*exp(v(V)/(1-v(V)))*ei(-(v(V)/(1-v(V)))-(NoB*v(V)/Si(k)))+NoB*ei(-(NoB*v(V)/Si(k)));
end
end

plot(K,10*log10(smooth(pav(:,2))),'-o');
hold on;
plot(K,10*log10(smooth(pavt(:,2))),'-*');
hold on;
plot(K,10*log10(smooth(pav(:,3))),'-o');
hold on
plot(K,10*log10(smooth(pavt(:,3))),'-*');
hold on
plot(K,10*log10(smooth(pav(:,4))),'-o');
hold on
plot(K,10*log10(smooth(pavt(:,4))),'-*');
xlabel('Qavg(dB)');
ylabel('Pavg(dB)');
title('Pavg vs Qavg - average received power constraint');


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

hold on;
plot(K,10*log10(Pavg(1,:)));
hold on;
plot(K,10*log10(Pavg(2,:)),'--')
hold on;
plot(K,10*log10(Pavg(3,:)),'-.')
axis([-5,5,-5,25])
legend({'var = 0.01 pav_er','var = 0.01 pav_er_th','var = 0.05 pav_er','var = 0.05 pav_er_th','var = 0.1 pav_er','var = 0.1 pav_er_th','var = 0.01 pav_tifr','var = 0.05 pav_tifr','var = 0.1 pav_tifr'},'Location','northwest');