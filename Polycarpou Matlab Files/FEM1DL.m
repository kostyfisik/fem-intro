%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NAME OF FILE: FEM1DL.m
%
% PURPOSE: This Matlab code solves a one-dimensional boundary-value problem
% using the finite element method. The boundary-value problem at hand is
% defined by the Poisson's equation in a simple medium (isotropic, linear,
% and homogeneous) characterized by a uniform electron charge density rhov=-rho0. 
% The finite element method is based on sub-dividing the one-dimensional domain 
% into Ne linear elements of equal length. Two Dirichlet boundary conditions are 
% applied: V=Va at the leftmost node of the domain, and V=Vb at the rightmost 
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
Ne=4; % Number of linear elements
% ==========================
eps=epsr*8.85*10^-12;
Nn=Ne+1;
for i=1:Ne
    elmconn(i,1)=i;
    elmconn(i,2)=i+1;
end
le=L/Ne;
Ke(1,1)=eps/le;
Ke(1,2)=-eps/le;
Ke(2,1)=-eps/le;
Ke(2,2)=eps/le;
fe=-le*rho0/2*([1;1]);
K=zeros(Nn);
f=zeros(Nn,1);
for e=1:Ne
    for i=1:2
        for j=1:2
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
    x(e,1)=(e-1)*le;
    x(e,2)=e*le;
end

Npoints=50;
dx=le/(Npoints);
for e=1:Ne
    for i=1:Npoints
        idx=(e-1)*Npoints+i;
        xeval(idx)=(idx-1)*dx;
        xeval2(idx+e-1)=(idx-1)*dx;
        Veval(idx)=V(e)*(x(e,2)-xeval(idx))/le+V(e+1)*(xeval(idx)-x(e,1))/le;
        Eeval(idx+e-1)=(V(e)-V(e+1))/le;
    end
    if e < Ne
        xeval2(idx+e)=xeval2(idx+e-1);
        Eeval(idx+e)=(V(e+1)-V(e+2))/le;
    end
end
xeval(idx+1)=idx*dx;
xeval2(idx+e)=idx*dx;
Veval(idx+1)=V(Ne+1);
Eeval(idx+e)=Eeval(idx+e-1);
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
ylabel('V (Volts)');
legend('FEM','Exact');
%
% Error Analysis (Using the area bounded between the two curves per element)
% =========================================================================
PercentError=0.0;
Aextot=0.0;
Diff=0.0;
for e=1:Ne
    Aexact=rho0/(6*eps)*(x(e,2)^3-x(e,1)^3)-0.5*(rho0*L/(2*eps)+Va/L)*(x(e,2)^2-x(e,1)^2)+Va*(x(e,2)-x(e,1));
    Aextot=Aextot+Aexact;
    Afe=le/2*(V(e)+V(e+1));
    Diff=Diff+abs(Aexact-Afe);
end
Aextot=abs(Aextot);
PercentError=Diff/Aextot*100
%
% Error Analysis (Using the L2-Norm Definition)
% ============================================
L2Error=0.0;
for e=1:Ne
    f1=rho0^2*(x(e,2)^5-x(e,1)^5)/(20*eps^2);
    f2=(V(e)/le-rho0*L/(2*eps)-Va/L-V(e+1)/le)*rho0/(4*eps)*(x(e,2)^4-x(e,1)^4);
    f3=(1/3)*((V(e+1)*x(e,1)/(eps*le)+Va/eps-V(e)*x(e,2)/(eps*le))*rho0+(V(e)/le-rho0*L/(2*eps)-Va/L-V(e+1)/le)^2)*(x(e,2)^3-x(e,1)^3);
    f4=(V(e+1)*x(e,1)/le+Va-V(e)*x(e,2)/le)*(V(e)/le-rho0*L/(2*eps)-Va/L-V(e+1)/le)*(x(e,2)^2-x(e,1)^2);
    f5=(V(e+1)*x(e,1)/le+Va-V(e)*x(e,2)/le)^2*(x(e,2)-x(e,1));
    L2Error=L2Error+(f1+f2+f3+f4+f5);
end
L2Error=sqrt(L2Error)


figure(2);
plot(xeval2,Eeval,'k--'); % Plot the Electric field obtained from the finite element method
hold
plot(xeval,Eexact,'k-'); % Plot the Electric field obtained from the exact analytical solution
xlabel('x (meters)');
ylabel('E (V/m)');
legend('FEM','Exact');


