function y = surfd(s,x)
r = 0.5;
y = 0.5 * s.alpha * (1.0 + x ./ (r^2 + x.^2)^(1.0/2.0));
end