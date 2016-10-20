function plotSin(y, y2)
%y is frequency
x=linspace(0,2*pi,y*16+1);
figure

if nargin==1
    plot(x, sin(y*x));
else
    disp('two inputs were given; what is wrong with you?')
end