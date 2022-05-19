clc; clear all; close all;

s.alpha =1.9;       % Scale angle
fprintf("Angle is %f \n",s.alpha)

s.I = 201;      % Number of streamwise cells
s.J = 101;      % Number of wall-normal cells
s.ymax = 50;    % Domain height

nt = 500;            % Iterations

s.a = 5;           % Streamwise stretching parameter
s.b = 5;           % Wall-normal stretching parameter

s = start(s);       % Initial conditions          

s1 = s;
    

tol1 = 1.0;  
tol = tol1;
iter = 0;
firstsec=1;
% while (tol > 5.0e-4)  % Convergence based on tolerance

% Record a movie
file = "/grid_convergence/movie_alpha" +num2str(s.alpha) + "I" + num2str(s.I) + "J" + num2str(s.J);
filename = strcat(pwd, file);
vidObj = VideoWriter(filename + '.avi');
open(vidObj)

while (iter < nt)       % Convergence based on number of iterations
    
s1 = update(s);         % Advance the solution
tol = norm(abs(s.tau(:,1) - s1.tau(:,1))) / tol1;    % L2 norm, normalised by the first value
if iter==0
    tol1 = tol;
end

s = s1;
iter = iter +1;
fprintf("Tol is %e at iteration %i \n", tol, iter)

% Record a movie

plot(s1.xr, s1.tau(:,1),'-b')
xlim([-20.0, 20.0])
ylim([-2, 1.5])
hold on
plot(s1.xr,zeros(s.I,1), 'r--')
hold off
xlabel('$x$','interpreter','latex','fontsize',15)
ylabel('$\tau$','interpreter','latex', 'fontsize', 15)

currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);

end
close(vidObj);
export(s1);            % Write surface data to .dat file, field data to .tec file

% figure(1)
plot(s1.xr, s1.tau(:,1))
xlim([-20.0, 20.0])
ylim([-1.0, 1.5])
hold on
plot(s1.xr,zeros(s.I,1), 'r--')
xlabel('$x$','interpreter','latex','fontsize',15)
ylabel('$\tau$','interpreter','latex', 'fontsize', 15)

