function write_tecplot(filename, s)

fileID = fopen(filename, 'w');

% Calculate shear stress
imax = s.I;
jmax = s.J;

fprintf(fileID, '"variables=\"x\",\"y\",\"u\",\"v\",\n"');
fprintf(fileID, "zone i=%i, j=%i, f=point\n", imax, jmax);

for i=1:imax
    for j=1:jmax
        fprintf(fileID, '%16.16f %16.16f %16.16f %16.16f \n',s.x(i), ...
            s.y(j), s.u(i,j), s.v(i,j));
    end
end

fclose(fileID);

end