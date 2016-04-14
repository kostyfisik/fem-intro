%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NAME OF FILE: FEM1DQ.m
%
% PURPOSE: This Matlab code solves a one-dimensional boundary-value problem
% using the finite element method. The boundary-value problem at hand is
% defined by the Poisson's equation in a simple medium (isotropic, linear,
% and homogeneous) characterized by a uniform electron charge density rhov=-rho0. 
% The finite element method is based on sub-dividing the one-dimensional domain 
% into Ne quadratic elements of equal length. Two Dirichlet boundary conditions 
% are applied: V=Va at the leftmost node of the domain, and V=Vb at the rightmost 
% node of the domain. The length of the domain is denoted by L. The user has 
% the freedom to modify any of the defined input variables including charge 
% density (rho0), domain length (L), boundary conditions (Va & Vb), dielectric 
% constant (epsr), and number of elements (Ne).
%
% All quantities are expressed in the SI system of units.
%
% The Primary Unknown Quantity is the Electric Potential.
% The Secondary Unknown Quantity is the Electric Field.
%
% Written by Anastasis Polycarpou (Last updated: 8/12/2005)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear % Clear all variables
%
% Define the Input Variables
% ==========================
L=8*10^-2; % Length of the domain
rho0=1*10^-8; % Charge density
epsr=1.0; % Dielectric constant of the domain
Va=1; % Boundary condition at the leftmost node
Vb=0; % Boundary condition at the rightmost node
Ne=4; % Number of quadratic elements
% ==========================
eps=epsr*8.85*10^-12;
Nn=2*Ne+1;
for e=1:Ne
    elmconn(e,1)=2*e-1;
    elmconn(e,2)=2*e+1;
    elmconn(e,3)=2*e;
end
le=L/Ne;
Ke(1,1)=7*eps/(3*le);
Ke(1,2)=eps/(3*le);
Ke(1,3)=-8*eps/(3*le);
Ke(2,1)=eps/(3*le);
Ke(2,2)=7*eps/(3*le);
Ke(2,3)=-8*eps/(3*le);
Ke(3,1)=-8*eps/(3*le);
Ke(3,2)=-8*eps/(3*le);
Ke(3,3)=16*eps/(3*le);
fe=-le*rho0/6*([1;1;4]);
K=zeros(Nn);
f=zeros(Nn,1);
for e=1:Ne
    for i=1:3
        for j=1:3
            K(elmconn(e,i),elmconn(e,j))=K(elmconn(e,i),elmconn(e,j))+Ke(i,j);
        end
        f(elmconn(e,i))=f(elmconn(e,i))+fe(i);
    end
end
for i=2:Nn
    f(i)=f(i)-K(i,1)*Va;
end
K(:,1)=0;
K(1,:)=0;
K(1,1)=1;
f(1)=Va;

for i=2:Nn-1
    f(i)=f(i)-K(i,Nn)*Vb;
end
K(:,Nn)=0;
K(Nn,:)=0;
K(Nn,Nn)=1;
f(Nn)=Vb;

V=K\f;
for e=1:Ne
    x(e,1)=(2*e-2)*le/2;
    x(e,2)=(2*e)*le/2;
    x(e,3)=(2*e-1)*le/2;
end

Npoints=10;
dx=le/(Npoints);
for e=1:Ne
    for i=1:Npoints
        idx=(e-1)*Npoints+i;
        xeval(idx)=(idx-1)*dx;
        ksi=2*(xeval(idx)-x(e,3))/le;
        Veval(idx)=V(2*e-1)*0.5*ksi*(ksi-1)+V(2*e+1)*0.5*ksi*(ksi+1)+V(2*e)*(1+ksi)*(1-ksi);
        Eeval(idx)=(-2/le)*(V(2*e-1)*(ksi-0.5)+V(2*e+1)*(ksi+0.5)+V(2*e)*(-2*ksi));
    end
end
xeval(idx+1)=idx*dx;
Veval(idx+1)=V(2*Ne+1);
ksi=2*(xeval(idx+1)-x(Ne,3))/le;
Eeval(idx+1)=(-2/le)*(V(2*Ne-1)*(ksi-0.5)+V(2*Ne+1)*(ksi+0.5)+V(2*Ne)*(-2*ksi));
plot(xeval,Veval,'k--'); % Plot the Electric potential obtained from the finite element method
%
% Exact Analytical Solution
% =========================
hold
for i=1:Ne*Npoints+1
    Vexact(i)=rho0/(2*eps)*xeval(i)^2+((Vb-Va)/L-rho0*L/(2*eps))*xeval(i)+Va;
    Eexact(i)=rho0*(L-2*xeval(i))/(2*eps)+(Va-Vb)/L;
end
plot(xeval,Vexact,'k-'); % Plot the Electric potential obtained from the exact analytical solution
xlabel('x (meters)');
ylabel('V (volts)');
legend('FEM','Exact');
%
% Error Analysis (Using the area bounded between the two curves per element)
% =========================================================================
PercentError=0.0;
Aextot=0.0;
Diff=0.0;
for e=1:Ne
    g1=rho0/(6*eps)*(x(e,3)^3-x(e,1)^3);
    g2=0.5*(rho0*L/(2*eps)+Va/L)*(x(e,3)^2-x(e,1)^2);
    g3=Va*(x(e,3)-x(e,1));
    A1exact=g1-g2+g3;
    A1fe=le/24*(5*V(2*e-1)-V(2*e+1)+8*V(2*e));
    Diff=Diff+abs(A1exact-A1fe);
    g1=rho0/(6*eps)*(x(e,2)^3-x(e,3)^3);
    g2=0.5*(rho0*L/(2*eps)+Va/L)*(x(e,2)^2-x(e,3)^2);
    g3=Va*(x(e,2)-x(e,3));
    A2exact=g1-g2+g3;
    A2fe=le/24*(-V(2*e-1)+5*V(2*e+1)+8*V(2*e));
    Diff=Diff+abs(A2exact-A2fe);
    Aextot=Aextot+A1exact+A2exact;
end
Aextot=abs(Aextot);
PercentError=Diff/Aextot*100

figure(2);
plot(xeval,Eeval,'k--'); % Plot the Electric field obtained from the finite element method
hold
plot(xeval,Eexact,'k-'); % Plot the Electric field obtained from the exact analytical solution
xlabel('x (meters)');
ylabel('E (V/m)');
legend('FEM','Exact');
