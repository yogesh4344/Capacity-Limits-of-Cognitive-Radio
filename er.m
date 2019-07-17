clc;
clear all;
close all;
NoB = 1;
for j=1:10^3
hs(j) = exprnd(1);
hp(j) = exprnd(1);
syms y;
P(j) = ((y/(hp(j)))-(NoB/hs(j)))*hp(j);
end
Q = mean(P) == 10^(-0.5);
S = solve(Q,y);
Si = double(S);
for i = 1:10^3
p = (Si/(hp(i)))-(NoB/hs(i));
C(i) = log(1+((p*hs(i))/NoB));
end
c = mean(C);
Ct = log((Si + NoB)/NoB);