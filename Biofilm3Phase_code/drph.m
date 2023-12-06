function f = drph(b,dr,m)
 f = 0*b;
%  i = 2:m;
%  f(1) = 0;
%  f(i) = (b(i)- b(i-1))./dr;

 i = 2:m-1;
 f(1) = 0;
 f(i) = (b(i+1) - b(i-1))./(2.*dr);
 f(m) = 0;
 
%  f = 0*b;
%  i = 3:m;
%  f(1) = 0;
%  f(2) = (b(2)- b(1))./dr;
%  f(i) =(3.*b(i) - 4.*b(i-1) + b(i-2))./(2*dr);
end