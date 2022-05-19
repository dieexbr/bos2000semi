function psi = find_psi(s)

dy = s.dy;
imax = s.I;
jmax = s.J;

for i =1:imax
    sum = 0.0;
    for j=2:jmax
        ga = omega(s.y(j));
        gb = omega(s.y(j-1));
        sum = sum + (dy/2.0) * (s.b * s.u(i,j)/ga + s.b * s.u(i,j-1)/gb);
        psi(i,j) = sum;
    end
end

end