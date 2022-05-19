function s = start(s);

s.xmin = -1; % xmin
s.xmax = 1;  % xmax

s.dx = (s.xmax - s.xmin) / (s.I-1);

% Initialise arrays;
s.x = linspace(s.xmin, s.xmax, s.I);    % Streamwise position in transformed coordinates
s.xr = s.a* tan(s.x * pi/2.0);          % Streamwise position in real coordinates

s.hatymax = 2.0 / pi * atan(s.ymax / s.b);  % Domain height in transformed coordinates

s.dy = (s.hatymax - 0) / (s.J-1);
s.y = linspace(0, s.hatymax, s.J);      % Wall normal position in transformed coordinates
s.yr = s.b * tan(s.y * pi/2.0);         % Wall normal position in real coordinates

s.tau = ones(s.I, s.J);
s.v = zeros(s.I, s.J);

% Initialise solution
for i=1:s.I
   s.u(i,:) = s.yr;
end
s.p = zeros(s.I);

end