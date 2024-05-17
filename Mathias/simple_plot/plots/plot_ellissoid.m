function plot_ellissoid(ba, c)

a = c.MajorAxisLength/2; 
b = c.MinorAxisLength/2; 
Xc = c.Centroid(1); 
Yc = c.Centroid(2);
phi = deg2rad(-c.Orientation);
t = linspace(0,2*pi,50);
x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
figure, imagesc(ba(:,500:end)); 
set(gca,'YDir','normal')
hold on 
plot(x,y,'r','Linewidth',1)