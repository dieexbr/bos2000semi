function y = surfd2(s,x)
r = 0.5;
y = r^2 * s.alpha./(2.0*((r^2 + x.^2).^(3/2)));
end