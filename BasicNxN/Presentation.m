clear; close;

ratio =  [2 3];
n = 1e3;
t_i = 0;
t_e = 2*pi;


q = ratio(2)/ratio(1);
dt = (t_e-t_i)/n;
t = t_i:dt:(t_e*lcm(ratio(1),ratio(2)));

a = 1;
b = 1.5;
c = sqrt(abs(b^2-a^2));

x = a*cos(t);
y = b*sin(t)-c;



x2 = x.*cos(t) - y.*sin(t);
y2 = x.*sin(t) + y.*cos(t);

% plot(x,y,'-k')
plot(x2,y2,'-k')
axis equal