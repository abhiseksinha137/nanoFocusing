function [dn]=dnBYdw(em,ed, k0,R)

c=3e8;
w0=k0*c;
dw=w0/10000;

k=@(w) w/c;

n1=newtonRaphson2(em,ed, k(w0),R)
n2=newtonRaphson2(em,ed, k(w0+dw),R)

dn= (n2-n1)/dw;