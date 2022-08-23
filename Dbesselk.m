function [ret]=Dbesselk(v, z)
ret= -besselk(v-1,z) - v./z.*besselk(v,z);
