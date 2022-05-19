function export(s)

file = "/grid_convergence/alpha" +num2str(s.alpha) + "I" + num2str(s.I) + "J" + num2str(s.J);

% Write wall data to file
filename = strcat(pwd, file);
[status, msg, msgID] = mkdir(filename);
filename = strcat(filename, "/results_surf.dat");

write_wall(filename, s);

% Write tecplot data to file
filename = strcat(pwd, file) ;
filename = strcat(filename, "/results.tec");
write_tecplot(filename, s);

% Save structure to restart next angle
filename = strcat("alpha", num2str(s.alpha) ,".mat");
save(filename, '-struct', 's');

end