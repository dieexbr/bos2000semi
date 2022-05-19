function write_wall(filename, s)

fileID = fopen(filename, 'w');

% Calculate shear stress
imax = s.I;
jmax = s.J;

for i=1:imax
    fprintf(fileID,"%16.16f %16.16f %16.16f \n", s.a*tan(s.x(i)*pi/2.0), s.tau(i,1), s.p(i));
end

fclose(fileID);

end