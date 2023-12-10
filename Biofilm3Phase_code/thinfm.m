function f = thinfm(t,c,dr,n,psin,psim,Q,gmma,kb,D, kmb, kms, pe, hstar)
% f = 0*c;
psid = 0; 

s0 = 6.66446; beta = 320; k1 = 0.5;

h = c(1:n); phic = c(n+1:2*n); phie = c(2*n+1:3*n); 
phia = c(3*n+1:4*n); cs = c(4*n+1:5*n); cb = c(5*n+1:6*n);
phim = 1 - phia; % phim = phic + phie;
hv = @(i) heaviside(h(i)- hstar); r = @(i) dr.*(i-1);

psix = 3.*(phie.^2).*drph(phie,dr,n);
th = tht(h,dr,n); %theta
dgb = drgb(cb,h,dr,n); dgs = drgs(cs,dr,n); 
phie0 = 0.04; v0 = Q.*(phie.^3 - phie0.^3);

uh = - (s0./3).*thinf(h.^3,1./phim, psix, dr, n) + (gmma./3).*thinf(h.^3, 1./phim, th, dr, n) ...
    + (gmma./beta).*thinf(h,(phia.^2)./phim, th, dr, n);
uc = - (s0./3).*thinf(h.^3, phic./phim, psix, dr, n) + (gmma./3).*thinf(h.^3, phic./phim, th, dr, n);
ue = - (s0./3).*thinf(h.^3, phie./phim, psix, dr, n) + (gmma./3).*thinf(h.^3, phie./phim, th, dr, n);
ua = - (s0./3).*thinf(h.^3, phia./phim, psix, dr, n) + (gmma./3).*thinf(h.^3, phia./phim, th, dr, n) ...
    + (gmma./beta).*thinf(h, (phia.^2)./phim, th, dr, n);
ucb = - (s0./3).*thinf(h.^3, cb./phim, psix, dr, n) + (gmma./3).*thinf(h.^3, cb./phim, th, dr, n) ...
    + (gmma./beta).*thinf(h,(cb.*phia.^2)./phim, th, dr, n);

gc = psin.*phic.*cb./(k1 + cb)  - psid.*phic; ge = psim.*phic.*cb./(k1 + cb)  + psid.*phic;
gh = (ge + gc).*h + phia.*v0;
for i = 1:n
f1 = gh(i); %rxn h
f2 = gc(i).*h(i);% + phic(i).*v0(i); %rxn phic
f3 = ge(i).*h(i);% + phie(i).*v0(i); %rxn phie
f4 = phia(i).*v0(i); %rxn phia

f5 = - hv(i).*D.*kms.*(cs(i) - cb(i)); %rxn gs
f6 = - kmb*(cb(i) - cs(i)) - kb.*phic(i).*cb(i).*h(i)./(k1 + cb(i)); %rxn gb

dht(i) =  - uh(i) + hv(i).*f1; dhdt = -uh(i) + f1;
dphict(i) = (- phic(i).*dhdt - uc(i) + f2)./h(i);
dphiet(i) = (- phie(i).*dhdt - ue(i) + f3)./h(i);
dphiat(i) = (- phia(i).*dhdt - ua(i) + f4)./h(i);
dgst(i) = D.*dgs(i) + f5;
dgbt(i) = hv(i).*(dgb(i) + f6 -  pe.*ucb(i))./(pe.*h(i));
end
f = [dht, dphict, dphiet, dphiat, dgst, dgbt]; f = f';
end