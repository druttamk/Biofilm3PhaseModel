function f = tht(h,dr,n)%theta = d((1/r)*d(r*dh/dr)/dr)/dr
f = 0.*h;
r = @(i) dr.*(i-1);
hf = @(b,i) (b(i+1) + b(i))./2;

for j =1
f(1) = 0;
end
for j = 2
f(2) = (r(3).*(h(4)-h(2))- r(2).*(h(3)-h(1)))/((r(2) + r(3)).*dr^3) ...
    - r(2).*(h(3)-h(1))/(r(2).*dr^3);
end

for j = 3:n-2
h_1 = r(j+1).*(hf(h,j+1) - hf(h,j));
h_2 = r(j).*(hf(h,j) - hf(h,j-1));
h_3 = r(j-1).*(hf(h,j-1) - hf(h,j-2));

h_f = (h_1 - h_2)./hf(r,j);
h_b = (h_2 - h_3)./hf(r,j-1);

f(j) = (h_f - h_b)./dr^3;
end

for j = n-1
f(n-1) = (r(n).*(3*h(n) - 4*h(n-1) + h(n-2)) - r(n-1).*(h(n)-h(n-2)))/((r(n) + r(n-1)).*dr^3) ...
    - (r(n-1).*(h(n)-h(n-2))- r(n-2).*(h(n-1)-h(n-3)))/((r(n-1) + r(n-2)).*dr^3);
end

for j = n
h1 = 3*r(n).*(3*h(n)-4*h(n-1)+h(n-2)) - 4*(r(n) + r(n-1))*(h(n)- h(n-1))...
    + r(n-1)*(h(n) - h(n-2));
h2 = r(n).*(3*h(n)-4*h(n-1)+h(n-2)) - r(n-1).*(h(n)-h(n-2));
h3 = (r(n) + r(n-1))*(h(n) - h(n-1)) - (r(n-1) + r(n-2))*(h(n-1) - h(n-2));

c1 = 3/(r(n)*dr^2); c2 = 8/((r(n) + r(n-1))*dr^2); c3 = 1/(r(n-1)*dr^2);
f(n) = (c1*h1 - c2*h2 +  c3*h3)/(2*dr) ;

% f(n) = ((3*h(n-1) - 4*h(n) + h(n))*(2*n/(2*n-1))...
%     - r(n-1)*(h(n) - h(n-2))/(2*dr*hf(r,n-1)))/dr^3;
end
end