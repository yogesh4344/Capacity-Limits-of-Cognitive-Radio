clc;
clear all;
close all;
NoB = 1;
v = [0,0.01,0.1,0.3];
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

for i = 1:2000
    hpc = exprnd(1-v(V));
    if hpc >= 0
        hpc = hpc;
    else hpc = 0;
    end
    hsc = exprnd(1);
    
    p(i) = (Si(k)/(hpc+v(V))) - NoB/hsc;
    if p(i) >= 0
        p(i) = p(i);
    else p(i) = 0; end
    C(i) = log(1+(p(i)*hsc/NoB));
    
end
c(k,V) = mean(C);
if v(V) == 0
    Ct(k,V) = log((Si(k)+NoB)/NoB);
else
Ct(k,V) = exp(v(V)/(1-v(V)))*ei(-(v(V)/(1-v(V)))-(NoB*v(V)/Si(k)))-ei(-(NoB*v(V)/Si(k)));
end
end
end
plot(K,smooth(c(:,1)),'-o');
hold on;
plot(K,Ct(:,1),'-*');
hold on;
plot(K,smooth(c(:,2)),'-o');
hold on
plot(K,Ct(:,2),'-*');
hold on
plot(K,smooth(c(:,3)),'-o');
hold on
plot(K,Ct(:,3),'-*');
hold on
plot(K,smooth(c(:,4)),'-o');
hold on
plot(K,Ct(:,4),'-*');
legend({'var = 0 cer','var = 0 cer_th','var = 0.01 cer','var = 0.01 cer_th','var = 0.1 cer','var = 0.1 cer_th','var = 0.3 cer','var = 0.3 cer_th'},'Location','northwest');
xlabel('Qavg (dB)');
ylabel('Cerg/B');
title('Cerg/B vs Qavg - average received power constraint');
