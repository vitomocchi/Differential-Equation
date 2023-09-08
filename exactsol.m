
% Finding the exact solution with syms tool.
syms u(x) a(x) f(x)
Du=diff(u,1);
f=1;
a=exp(-x);
temp_uex=dsolve(-diff((a*Du),1)==f);
temp_uex=simplify(temp_uex); 
% uex=vectorize(uex)
temp_uex
temp_duex=diff(temp_uex,1)
der_value_in0=subs(temp_duex,x,0)
C1=0 % value of integration constant. Can be seen as u'(0).
uex = exp(x)*(C1 - x + 1);
duex=diff(uex,1)
der_value_in1=subs(duex,x,1)