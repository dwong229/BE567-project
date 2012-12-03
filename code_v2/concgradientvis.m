

% Concentration gradient visualization

% Exponential decay of attractant

% Attractant Source at [xs,ys]
clear
%xs = -20e4;
%ys = 20e4;

xs = 300;
ys = 300;

%xcoord = linspace(-21e4,5e4,10);
%ycoord = linspace(-5e4,21e4,10);

xcoord = -500:50:500;
ycoord = -500:50:500;


for x = 1:length(xcoord)
    concgrad(:,x) = 100*exp(-sqrt((xs-xcoord(x))^2+(ys-ycoord).^2)/300);
end
%contourf(xcoord,ycoord,concgrad)
pcolor(xcoord,ycoord,concgrad)
%mesh(xcoord,ycoord,concgrad)        

%concgrad = arrayfun(@(xcoord,ycoord)(exp(-sqrt((xs-xcoord)^2+(ys-ycoord)^2))),xcoord,ycoord);
shading flat
shading interp
