function v = find_v(s);

% v on the bottom
imax = s.I;
jmax = s.J;
dx = s.dx;
dy = s.dy;

psi = s.psi;

v = zeros(imax,jmax);

v(1,:) = zeros(1,jmax);
v(imax,:) = zeros(1,jmax);
v(:,1) = zeros(imax,1);

val1 = - omega(s.x(1)) / s.a;
valI = - omega(s.x(imax)) / s.a;

for j=2:jmax
    v(1,j) = val1 * (-3.0 * psi(1,j) + 4.0 * psi(2,j) - psi(3,j))/(2.0 * dx);
    v(imax,j) = valI * (3.0 * psi(imax,j) - 4.0 * psi(imax-1,j) + psi(imax-2,j))/(2.0*dx);
end

for i=2:imax-1
    for j=2:jmax
        v(i,j) = - omega(s.x(i))/s.a * (psi(i+1,j) - psi(i-1,j)) / (2.0 * dx);
    end
end

end