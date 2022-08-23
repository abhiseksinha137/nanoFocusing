function [ret]=K(n,x)
ret=pi/2 * (I(-n,x) - I(n,x))/sin(n*pi);