function y = surf(s, x)
r = 0.5;
y = 0.5 * s.alpha .* (x + sqrt(x.^2 + r.^2));
end