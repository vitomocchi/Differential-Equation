clearvars, close all

% % Finding the exact solution with syms tool.
% syms u(x) a(x) f(x)
% Du=diff(u,1);
% f=1;
% a=exp(-x);
% temp_uex=dsolve(-diff((a*Du),1)==f);
% temp_uex=simplify(temp_uex); 
% %uex=vectorize(uex)
% temp_uex
% temp_duex=diff(temp_uex,1)
% der_value_in0=subs(temp_duex,x,0)
% C1=0 % value of integration constant. Can be seen as u'(0).
% uex = exp(x)*(C1 - x + 1);
% duex=diff(uex,1)
% der_value_in1=subs(duex,x,1)

L=1; % b-a
n=10; % # spacings 
h=L/(n+1); % mesh step (n+1) elements, (n+2) total nodes.
f=@(x)1+0.*x; 
a=@(x)exp(-x);

g0=0; g1=-exp(1); % Neumann BCs.

uex=@(x)(exp(x).*(g0 - x + 1));

xnodes=linspace(0,L,n+2);

diag = [-3/2, a(xnodes(2:end-1)-1/2*h)+a(xnodes(2:end-1)+1/2*h) 3/2]; % Added first and last element to implement one-sided derivatives approximation.
offdiaginf = [-a(xnodes(2:end-1)-1/2*h), -2, 0]; % spdiags ignores last value.
offdiagsup = [0, 2, -a(xnodes(2:end-1)+1/2*h)]; % spdiags ignores first value.

% Create the sparse matrix using spdiags().
A = spdiags([offdiaginf', diag', offdiagsup'], -1:1, n+2, n+2);

% The following values need to be add in order tu fully implement Neumann conditions (using the one-sided derivatives approximation).
A(1,3)=-1/2;
A(end,end-2)=1/2; 

A=A/h^2;

b=[g0/h, f(xnodes(2:end-1)), g1/h]; b=b'; % known term (column) + Neumann conditions.

Uh=A\b;

subplot(2,1,1);
plot(xnodes,Uh,'b--')
xlabel('x'), ylabel('Uh(x)'), title('Approximated solution')
subplot(2,1,2);
plot(xnodes,uex(xnodes),'k-')
xlabel('x'), ylabel('uex(x)'), title('Exact solution')
fprintf('h = %12.10e err = %12.10e\n',h,norm(Uh'-uex(xnodes),'inf'));
