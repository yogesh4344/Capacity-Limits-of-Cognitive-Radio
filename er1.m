clc;
clear all;
close all;
NoB = 1;
v = 0.01;
K = -5:1:5;
for k = 1:length(K)
for j=1:10^3
hs(j) = exprnd(1);
hp(j) = exprnd(1-v);

syms y;
P(j) = ((y/(hp(j) + v))-(NoB/hs(j)))*hp(j) + ((y/(hp(j)))-(NoB/hs(j)))*v;
end
Q = mean(P) == 10^(K(k)/10);
S = solve(Q,y);
Si = double(S);
for i = 1:10^3
p = (Si/(hp(i) + v))-(NoB/hs(i));
if p >= 0
    p = p;
else
    p = 0;
end
C(i) = log(1+((p*hs(i))/NoB));
end
c(k) = mean(C);
Ct(k) = exp(-v/(1-v))*ei(-(v/(1-v)))-ei(-(NoB*v/Si));
end