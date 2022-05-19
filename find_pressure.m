function p = find_pressure(s)

tau = s.tau;

dx = s.dx;
dy = s.dy;

imax = s.I;
jmax = s.J;

p = [];

p(1) = 0.0;
p(imax) = s.alpha;

for i=2:imax-1
    omx = omega(s.x(i))/s.a;
    summin = 0;
    summax = 0;
    for j=2:jmax
       summin = summin + (tau(i-1,j) * s.b / omega(s.y(j)) -1) + (tau(i-1,j-1) * s.b / omega(s.y(j-1)) -1);
       summax = summax + (tau(i+1,j) * s.b / omega(s.y(j)) -1) + (tau(i+1,j-1) * s.b / omega(s.y(j-1)) -1);
    end
    summin = (dy/2.0) / (2.0 * dx) * summin;
    summax = (dy/2.0) / (2.0 * dx) * summax;
    
   p(i) = surfd(s,s.xr(i))  - omx * (summax - summin);
end

end
