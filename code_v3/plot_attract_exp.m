function [] = plot_attract_exp(fighandle,path)

figure = fighandle;
% figure out x/y axis bounds:
hold on

[~,xs,ys] = attract_exp([0,0]);
xmin = path(1)-10;
xmax = path(2)+10;
ymin = path(3)-10;
ymax = path(4)+10;
%xcoord = linspace(min(xs,path(1))-10,max(xs,path(2))+10,50);
%ycoord = linspace(min(ys,path(3))-10,max(ys,path(4))+10,50);
xcoord = linspace(path(1)-10,path(2)+10,50);
ycoord = linspace(path(3)-10,path(4)+10,50);
% copied from attract_exp.m
for x = 1:length(xcoord)
    concgrad(:,x) = 100*exp(-sqrt((xs-xcoord(x))^2+(ys-ycoord).^2)/300);
end
%contourf(xcoord,ycoord,concgrad)
pcolor(xcoord,ycoord,concgrad)
axis([xmin, xmax,ymin,ymax])
%mesh(xcoord,ycoord,concgrad)        

%concgrad = arrayfun(@(xcoord,ycoord)(exp(-sqrt((xs-xcoord)^2+(ys-ycoord)^2))),xcoord,ycoord);
shading flat
shading interp
