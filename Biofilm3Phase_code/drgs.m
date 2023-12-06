function f = drgs(gs,dr,m) %(1/r)*d(r*dgs/dr)/dr
f = 0*gs;
r = @(i) dr.*(i-1);
hf = @(b,i) (b(i+1) + b(i))./2;

for j = 1
f(j) = 2.*(gs(j+1) - gs(j))/(dr^2);
end

for j = 2:m-1
    rm = r(j).*(dr.^2);
f(j) = (hf(r,j).*(gs(j+1)- gs(j)) ...
    - hf(r,j-1).*(gs(j)-gs(j-1)))./rm;
end

for j = m
f(j) = 2*(gs(j-1) - gs(j))/(dr^2);
end
end