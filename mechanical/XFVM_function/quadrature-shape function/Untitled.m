%{
coord=[-1,-1;1,-1;1,1;-1,1];

xi=coord(1); eta=coord(2);
N=1/4*[ (1-xi)*(1-eta);
        (1+xi)*(1-eta);
        (1+xi)*(1+eta);
        (1-xi)*(1+eta)];
    %}

r=linspace(0,2,50);
theta=linspace(-180,180,50);
a=r'*cosd(theta);
b=r'*sind(theta);
%[x,y]=ndgrid(a,b);
[x1,y1]=ndgrid(r,theta);
z=sqrt(x1).*sind(y1/2);
surf(a,b,z)
%plot3(a,b,z)