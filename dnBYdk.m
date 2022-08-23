function [dn]=dnBYdk(em,ed, k0,R)

c=3e8;
dk=1e-3;


n1=newtonRaphson2(em,ed, k0,R);
n2=newtonRaphson2(em,ed, k0+dk,R);

dn= (n2-n1)/dk;