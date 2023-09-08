clearvars, close all

% Finding the exact solution with syms tool.
% syms u(x) a(x) f(x)
% Du=diff(u,1);
% f=1;
% a=exp(-x);
% uex=dsolve(-diff((a*Du),1)==f,u(0)==0,u(1)==0);
% uex=simplify(uex); 
% uex=vectorize(uex)

L=1; % b-a
n=10; % # spacings 
h=L/(n+1); % mesh step (n+1) elements, (n+2) total nodes.
f=@(x)1+0.*x; 
a=@(x)exp(-x);

uex=@(x)exp(x).*(1./(exp(1) - 1) - x + 1) - exp(1)./(exp(1) - 1);

xnodes=linspace(0,L,n+2);

u0=0; uL=0; % Dirichlet BCs.

diag = a(xnodes(2:end-1)-1/2*h)+a(xnodes(2:end-1)+1/2*h);
offdiaginf = [-a(xnodes(3:end-1)-1/2*h), 0]; % spdiags ignores last value.
offdiagsup = [0, -a(xnodes(2:end-2)+1/2*h)]; % spdiags ignores first value.

% Create the sparse matrix using spdiags().
A = spdiags([offdiaginf', diag', offdiagsup'],-1:1, n, n);

A=A/h^2;

b=f(xnodes(2:end-1)); b=b'; % known term (column).
b(1)=b(1)+u0*a(xnodes(2)-1/2*h)/h^2;  % modifies with Dirichlet BC in 0.
b(end)=b(end)+uL*a(xnodes(end-1)+1/2*h)/h^2; % modifies with Dirichlet BC in 1.

Uh=A\b; % dim Uh: (nx1)
Uh=[u0;Uh;uL];
plot(xnodes,Uh,'b*',xnodes,uex(xnodes), 'k--')
xlabel('x'), title('Exact solution and approximated solution: a comparison')
legend('Uh(x)','uex(x)')
fprintf('h = %12.10e err = %12.10e\n',h,norm(Uh'-uex(xnodes),'inf'));
