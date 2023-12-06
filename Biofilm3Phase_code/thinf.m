function f = thinf(h, phi, th, dr, m)%(1/r)*d(r*h*phi*theta)/dr
f = 0*h;
r = @(i) dr.*(i-1); hf = @(b,i) (b(i+1) + b(i))./2;

for j = 1
    f(j) = (4.*h(2).*phi(2).*th(2) - h(3).*phi(3).*th(3))./(2*dr);
end
for j = 2:m-1
    rj = r(j).*dr;
    f(j) = (hf(r,j).*hf(phi,j).*hf(h,j).*hf(th,j)...
        - hf(r,j-1).*hf(phi,j-1).*hf(h,j-1).*hf(th,j-1))./rj;
end
for j = m
    rm = r(j).*(2*dr);
    f(j) = (3*r(j).*h(j).*phi(j).*th(j)...
        - 4*r(j-1).*h(j-1).*phi(j-1).*th(j-1)...
        + r(j-2).*h(j-2).*phi(j-2).*th(j-2))./rm;
    
% %     f(j) = (- 4.*phi(j-1).*h(j-1).*th(j-1)...
% %         + phi(j-2).*h(j-2).*th(j-2))./(2.*dr);
end
end