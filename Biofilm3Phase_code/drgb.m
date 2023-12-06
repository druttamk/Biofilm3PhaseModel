function f = drgb(gb,h,dr,m) %(1/r)*d(r*h*dgb/dr)/dr
f = 0*gb;
r = @(i) dr.*(i-1);
hf = @(b,i) (b(i+1) + b(i))./2;

for j = 1
f(j) = 2.*hf(h,j).*(gb(2)- gb(1))/(dr.^2);
end

for j = 2:m-1
    r0 = r(j).*(dr.^2);
f(j) = (hf(r,j).*hf(h,j).*(gb(j+1)- gb(j))...
    - hf(r,j-1).*hf(h,j-1).*(gb(j)- gb(j-1)))./r0;
end

for j = m
f(j) = 2.*hf(h,j-1).*(gb(j-1) - gb(j))./dr^2;
end
end