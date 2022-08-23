% Heaviside Function
function[y]= THETA(x)

debug=false;
if (nargin < 1)
    debug=true;
    x=linspace(-10,10,1000);
end


y=ones(size(x));
y(x<0)=0;

if debug
    close all;
    plot(x,y, 'k','linewidth',2); title('Heaviside Function');
end