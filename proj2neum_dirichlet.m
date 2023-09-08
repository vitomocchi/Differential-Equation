clearvars, close all

% % Finding the exact solution with syms tool
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

% We now want to force a Dirichlet BC in addition to the Neumann BCs, in
% order to have an univocal solution.

L=1; % b-a
n=10; % # spacings 
h=L/(n+1); % Mesh step, (n+1) elements, (n+2) total nodes.
f=@(x)1+0.*x; 
a=@(x)exp(-x);

xnodes=linspace(0,L,n+2);

g0=0; g1=-exp(1); % Neumann BCs.
u0=1; % Initial condition, needed to uniform the approximated solution and the exact one - otherwise I'd have infinite solutions respecting Neumann BCs, different by a constant C. 

uex=@(x)(exp(x).*(g0 - x + 1));

diag = [a(xnodes(2:end-1)-1/2*h)+a(xnodes(2:end-1)+1/2*h), 3/2]; % Added first and last element to implement one-sided derivatives approximation.
offdiaginf = [-a(xnodes(3:end-1)-1/2*h), -2, 0]; % spdiags ignores last value.
offdiagsup = [0, -a(xnodes(2:end-1)+1/2*h)]; % spdiags ignores first value.

% Create the sparse matrix using spdiags().
A = spdiags([offdiaginf', diag', offdiagsup'], -1:1, n+1, n+1);

% The following value needs to be add in order tu fully implement Neumann conditions (using the one-sided derivatives approximation).
A(end,end-2)=1/2; 

A=[A; 2, -1/2, zeros(1,n-1)]; % Adding Neumann conditions in 0 (last row).

A=A/h^2;

b=[f(xnodes(2:end-1)), g1/h, g0/h+3/2*u0/h^2]; % Known term with Neumann conditions in x=0 (last position) and x=1 (penultimate position).
b(1)=b(1)+u0*a(xnodes(2)-1/2*h)/h^2;  % Modified with Dirichlet BC in 0 in order to have an unique solution.
b=b'; % Known term (column)

Uh=A\b;

Uh=[u0;Uh]; % Adding initial condition u0 for plotting.

subplot(2,1,1);
plot(xnodes,Uh,'b-')
xlim([0 1])
ylim([0 1])
xlabel('x'), ylabel('Uh(x)'), title('Approximated solution')
subplot(2,1,2);
plot(xnodes,uex(xnodes),'k-')
xlabel('x'), ylabel('uex(x)'), title('Exact solution')
fprintf('h = %12.10e err = %12.10e\n',h,norm(Uh'-uex(xnodes),'inf'));
