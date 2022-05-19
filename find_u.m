function u = find_u(s);
% Dependent on tau

imax = s.I;
jmax = s.J;
dx = s.dx;
dy = s.dy;

tau = s.tau;

u = zeros(imax,jmax);

u(1,:) = s.yr;
u(imax,:) = s.yr;
u(:,1) = zeros(imax,1);

for i=2:imax-1
    sum = 0.0;
    for j=2:jmax
        ga = omega(s.y(j));
        gb = omega(s.y(j-1));
        sum = sum + s.b *  (dy / 2.0) * ( s.tau(i,j) / ga + s.tau(i,j-1) / gb);
        u(i,j) = sum;
    end
end

end