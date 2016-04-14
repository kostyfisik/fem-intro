%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NAME OF FILE: FEM2DL_Cyl.m
%
% PURPOSE: This Matlab code solves a two-dimensional scattering problem
% using the finite element method. The scatterer is a perfectly
% conducting circular cylinder of radius robj; all the dimensions are given
% in terms of the free-space wavelength. A TM-to-z incident plane wave is 
% scattered from the circular cylinder and propagates undisturbed toward 
% infinity. To simulate this undisturbed propagation of the scattered field
% toward infinity, a first-order absorbing boundary condition (ABC) was
% imposed at a distance rho away from the center of the cylinder. The farther
% the ABC boundary is placed, the more accurate the ABC is. The free space 
% between the circular cylinder and the ABC boundary is subdivided into 
% triangles governed by linear interpolation functions. Dirichlet boundary 
% conditions are applied on the surface of the cylinder; i.e., the tangential 
% electric field, in this case Ez, is set to zero for all the nodes that 
% coincide with the surface of the circular cylinder.
%
% The finite element method is applied to the homogeneous scalar wave
% equation, otherwise known as the homogeneous Helmholtz equation. The
% primary unknown quantity is the total electric field in the z-direction
% which is given by the incident field plus the scattered field. The
% direction of the incident field is set toward the positive x-axis (phi_i=0) 
% whereas the total field is evaluated at a distance half-way between the 
% scatterer and the ABC boundary for all observation angles between 0 and 360 
% degrees. The numerical solution is compared with the exact analytical 
% solution.
%
% The user is allowed to set the following input parameters:
%
% rhoj = radius of the scatterer (circular cylinder) in wavelengths
% rho  = radius of the ABC boundary in wavelengths
% h    = discritization size in wavelengths
% E0   = amplitude of the incident electric field
%
% IMPORTANT: Depending on the number of nodes in the domain and the clock 
% speed of your computer, the finite element code may take several minutes 
% to execute. Try not to exceed 5,000 nodes otherwise you have to wait a 
% significant amount of time to get results!
%
% Written by Anastasis Polycarpou (Last updated: 8/12/2005)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear % Clear all variables
%
% Define the Input parameters 
% ===========================
robj=0.5; % Set the radius of the circular PEC cylinder in wavelengths
rho=1.5; % Set the radius of the outer circular ABC boundary in wavelengths
h=0.04; % Set the discretization size in wavelengths
E0=1; % Set the amplitude of the Incident field
% ===========================
%
% Determine the number of annular rings
% =====================================
Nseg=ceil((rho-robj)/(sqrt(3)*h/2));
dr=(rho-robj)/Nseg;
%
% Create the nodes on each annular ring
% =====================================
Nnods=0;
for J=1:Nseg+1
    r=robj+(J-1)*dr;
    Npoints=ceil(2*pi*r/h);
    dphi=2*pi/Npoints;
    for JJ=1:Npoints
        Nnods=Nnods+1;
        phi=(JJ-1)*dphi;
        x(Nnods)=r*cos(phi);
        y(Nnods)=r*sin(phi);
    end
end
%
% Triangulate using Delaunay triangulation
% ========================================
tri=delaunay(x,y);
%
% Eliminate the triangles within the perfectly conducting circular cylinder
% =========================================================================
Nelms=0;
for J=1:length(tri)
    Innods=0;
    for JJ=1:3
        r0=sqrt(x(tri(J,JJ))^2+y(tri(J,JJ))^2);
        if abs(robj-r0) < 0.0001*robj
            Innods=Innods+1;
        end
    end
    if Innods <= 2
        Nelms=Nelms+1;
        elmnod(Nelms,:)=tri(J,:);
    end
end
%
% Print the number of nodes and the number of elements
% ====================================================
fprintf('The number of nodes is: %6i\n',Nnods)
fprintf('The number of elements is: %6i\n',Nelms)
%
% Set all triangles in a counter-clockwise sense of numbering
% ===========================================================
for e=1:Nelms
    x21=x(elmnod(e,2))-x(elmnod(e,1));
    x31=x(elmnod(e,3))-x(elmnod(e,1));
    y21=y(elmnod(e,2))-y(elmnod(e,1));
    y31=y(elmnod(e,3))-y(elmnod(e,1));
    Ae=0.5*(x21*y31-x31*y21);
    if Ae < 0
        temp=elmnod(e,2);
        elmnod(e,2)=elmnod(e,3);
        elmnod(e,3)=temp;
    end
end
%
% Plot the mesh
% =============
figure(1)
triplot(elmnod,x,y);
xlabel('x (wavelengths)');
ylabel('y (wavelengths)');
axis([-rho rho -rho rho]);
axis square;
%
% Set the Dirichlet & Absorbing BCs on the cylinder and outer boundary,
% respectively
% =====================================================================
NEBC=0;
NABC=0;
for JJ=1:Nnods
    r0=sqrt(x(JJ)^2+y(JJ)^2);
    if abs(robj-r0) < 0.0001*robj
        NEBC=NEBC+1;
        EBC(NEBC)=JJ;
    elseif abs(rho-r0) < 0.0001*rho
        NABC=NABC+1;
        ABC(NABC)=JJ;
    end
end

for e=1:Nelms
    for I=1:3
        Nd=elmnod(e,I);
        flag(e,I)=0;
        for J=1:NABC
            if Nd == ABC(J)
                flag(e,I)=1;
            end
        end
    end
end 
%
% Print the number of Dirichlet BCs and the number of ABCs
% ========================================================
fprintf('The number of nodes on a Dirichlet boundary is: %6i\n',NEBC)
fprintf('The number of nodes on an ABC boundary is: %6i\n',NABC)
%
% Assign the constants for each element
% =====================================
mur0=1;
epsr0=1;
k0=2*pi;
gamma=k0*j+1/(2*rho);
for e=1:Nelms
    alphax(e)=1/mur0;
    alphay(e)=1/mur0;
    beta(e)=k0^2*epsr0;
    g(e)=0;
end
%
% Initialization of the global K matrix and right-hand side vector
% ================================================================
K=sparse(Nnods,Nnods);
b=zeros(Nnods,1);
%
% Form the element matrices and assemble to the global matrix
% ===========================================================
for e=1:Nelms
    x21=x(elmnod(e,2))-x(elmnod(e,1));
    x31=x(elmnod(e,3))-x(elmnod(e,1));
    x32=x(elmnod(e,3))-x(elmnod(e,2));
    x13=x(elmnod(e,1))-x(elmnod(e,3));
    y12=y(elmnod(e,1))-y(elmnod(e,2));
    y21=y(elmnod(e,2))-y(elmnod(e,1));
    y31=y(elmnod(e,3))-y(elmnod(e,1));
    y23=y(elmnod(e,2))-y(elmnod(e,3));
    Ae=0.5*(x21*y31-x31*y21);
    
    % Evaluation of the element K matrix
    % ==================================
    Me(1,1)=-(alphax(e)*y23^2+alphay(e)*x32^2)/(4*Ae);
    Me(1,2)=-(alphax(e)*y23*y31+alphay(e)*x32*x13)/(4*Ae);
    Me(2,1)=Me(1,2);
    Me(1,3)=-(alphax(e)*y23*y12+alphay(e)*x32*x21)/(4*Ae);
    Me(3,1)=Me(1,3);
    Me(2,2)=-(alphax(e)*y31^2+alphay(e)*x13^2)/(4*Ae);
    Me(2,3)=-(alphax(e)*y31*y12+alphay(e)*x13*x21)/(4*Ae);
    Me(3,2)=Me(2,3);
    Me(3,3)=-(alphax(e)*y12^2+alphay(e)*x21^2)/(4*Ae);
    
    % Evaluation of the element T matrix
    % ==================================
    Te(1,1)=beta(e)*Ae/6;
    Te(1,2)=beta(e)*Ae/12;
    Te(2,1)=beta(e)*Ae/12;
    Te(1,3)=beta(e)*Ae/12;
    Te(3,1)=beta(e)*Ae/12;
    Te(2,2)=beta(e)*Ae/6;
    Te(2,3)=beta(e)*Ae/12;
    Te(3,2)=beta(e)*Ae/12;
    Te(3,3)=beta(e)*Ae/6;
    
    % Sum the element matrices Me and Te
    % ==================================
    for I=1:3
        for J=1:3
            Ke(I,J)=Me(I,J)+Te(I,J);
        end
    end
    
    % Evaluation of element vector ge
    % ===============================
    ge(1)=g(e)*Ae/3;
    ge(2)=g(e)*Ae/3;
    ge(3)=g(e)*Ae/3;
    
    % Evaluation of the element vector pe & update of the element K-matrix
    % ====================================================================
    pe=zeros(3,1);
    if flag(e,1) == 1 && flag(e,2) == 1
        x1=x(elmnod(e,1));
        y1=y(elmnod(e,1));
        x2=x(elmnod(e,2));
        y2=y(elmnod(e,2));
        x21=x2-x1;
        Ledg=sqrt((x2-x1)^2+(y2-y1)^2);
        q0=gamma-j*k0*(y2-y1)/Ledg;
        C0=E0*q0*Ledg*exp(-j*k0*x1);
        if(x21 ~= 0)
            C1=(1-j*k0*x21-exp(-j*k0*x21))/(k0^2*x21^2);
            C2=(-1+j*k0*x21*exp(-j*k0*x21)+exp(-j*k0*x21))/(k0^2*x21^2);
        else
            C1=0.5;
            C2=0.5;
        end
        pe(1)=C0*C1;
        pe(2)=C0*C2;
        pe(3)=0;
        
        Ke(1,1)=Ke(1,1)-gamma*Ledg/3;
        Ke(1,2)=Ke(1,2)-gamma*Ledg/6;
        Ke(2,1)=Ke(2,1)-gamma*Ledg/6;
        Ke(2,2)=Ke(2,2)-gamma*Ledg/3;
    elseif flag(e,1) == 1 && flag(e,3) == 1
        x1=x(elmnod(e,1));
        y1=y(elmnod(e,1));
        x3=x(elmnod(e,3));
        y3=y(elmnod(e,3));
        x13=x1-x3;
        Ledg=sqrt((x1-x3)^2+(y1-y3)^2);
        q0=gamma-j*k0*(y1-y3)/Ledg;
        C0=E0*q0*Ledg*exp(-j*k0*x3);
        if(x13 ~= 0)
            C1=(1-j*k0*x13-exp(-j*k0*x13))/(k0^2*x13^2);
            C2=(-1+j*k0*x13*exp(-j*k0*x13)+exp(-j*k0*x13))/(k0^2*x13^2);
        else
            C1=0.5;
            C2=0.5;
        end
        pe(1)=C0*C2;
        pe(2)=0;
        pe(3)=C0*C1;
        
        Ke(1,1)=Ke(1,1)-gamma*Ledg/3;
        Ke(1,3)=Ke(1,3)-gamma*Ledg/6;
        Ke(3,1)=Ke(3,1)-gamma*Ledg/6;
        Ke(3,3)=Ke(3,3)-gamma*Ledg/3;
    elseif flag(e,2) == 1 && flag(e,3) == 1
        x2=x(elmnod(e,2));
        y2=y(elmnod(e,2));
        x3=x(elmnod(e,3));
        y3=y(elmnod(e,3));
        x32=x3-x2;
        Ledg=sqrt((x3-x2)^2+(y3-y2)^2);
        q0=gamma-j*k0*(y3-y2)/Ledg;
        C0=E0*q0*Ledg*exp(-j*k0*x2);
        if(x32 ~= 0)
            C1=(1-j*k0*x32-exp(-j*k0*x32))/(k0^2*x32^2);
            C2=(-1+j*k0*x32*exp(-j*k0*x32)+exp(-j*k0*x32))/(k0^2*x32^2);
        else
            C1=0.5;
            C2=0.5;
        end
        pe(1)=0;
        pe(2)=C0*C1;
        pe(3)=C0*C2;
        
        Ke(2,2)=Ke(2,2)-gamma*Ledg/3;
        Ke(2,3)=Ke(2,3)-gamma*Ledg/6;
        Ke(3,2)=Ke(3,2)-gamma*Ledg/6;
        Ke(3,3)=Ke(3,3)-gamma*Ledg/3;
    end   
    
    % Assemble element matrices & vectors into the global K matrix and b
    % vector
    % ==================================================================
    for I=1:3
        for J=1:3
            K(elmnod(e,I),elmnod(e,J))=K(elmnod(e,I),elmnod(e,J))+Ke(I,J);
        end
        b(elmnod(e,I))=b(elmnod(e,I))+ge(I)-pe(I);
    end
end
%
% Imposition of Dirichlet boundary conditions
% ===========================================
for I=1:NEBC
    for J=1:Nnods
        if J ~= EBC(I) 
            b(J)=b(J)-K(J,EBC(I))*0;
        end
    end
    K(:,EBC(I))=0;
    K(EBC(I),:)=0;
    K(EBC(I),EBC(I))=1;
    b(EBC(I))=0;
end
%
% Solution of the global matrix system
% ====================================
Ez=K\b;
%
% Generate the solution over a grid and plot it
% =============================================
%[xgrid,ygrid]=meshgrid(-rho:0.0025*(2*rho):rho,-rho:0.0025*(2*rho):rho);
[xgrid,ygrid]=meshgrid(-rho:0.01*(2*rho):rho,-rho:0.01*(2*rho):rho);
Ezgrid=zeros(101,101); %Ezgrid=zeros(401,401);
for I=1:101 %401
    for J=1:101 %401
       % I
       % J
        for e=1:Nelms
         
            x2p=x(elmnod(e,2))-xgrid(I,J);
            x3p=x(elmnod(e,3))-xgrid(I,J);
            y2p=y(elmnod(e,2))-ygrid(I,J);
            y3p=y(elmnod(e,3))-ygrid(I,J);
            A1=0.5*abs(x2p*y3p-x3p*y2p);
                
            x2p=x(elmnod(e,2))-xgrid(I,J);
            x1p=x(elmnod(e,1))-xgrid(I,J);
            y2p=y(elmnod(e,2))-ygrid(I,J);
            y1p=y(elmnod(e,1))-ygrid(I,J);
            A2=0.5*abs(x2p*y1p-x1p*y2p);
                
            x1p=x(elmnod(e,1))-xgrid(I,J);
            x3p=x(elmnod(e,3))-xgrid(I,J);
            y1p=y(elmnod(e,1))-ygrid(I,J);
            y3p=y(elmnod(e,3))-ygrid(I,J);
            A3=0.5*abs(x1p*y3p-x3p*y1p);
                
            x21=x(elmnod(e,2))-x(elmnod(e,1));
            x31=x(elmnod(e,3))-x(elmnod(e,1));
            y21=y(elmnod(e,2))-y(elmnod(e,1));
            y31=y(elmnod(e,3))-y(elmnod(e,1));
            Ae=0.5*(x21*y31-x31*y21);
                
            if abs(Ae-(A1+A2+A3)) < 0.00001*Ae   
                 ksi=(y31*(xgrid(I,J)-x(elmnod(e,1)))-x31*(ygrid(I,J)-y(elmnod(e,1))))/(2*Ae);
                 ita=(-y21*(xgrid(I,J)-x(elmnod(e,1)))+x21*(ygrid(I,J)-y(elmnod(e,1))))/(2*Ae);
                 N1=1-ksi-ita;
                 N2=ksi;
                 N3=ita;
                 Ezgrid(I,J)=N1*Ez(elmnod(e,1))+N2*Ez(elmnod(e,2))+N3*Ez(elmnod(e,3));
            end
        end
    end
end
%
% Plot the total electric field obtained from the finite element solution
% on a contour plot
% =======================================================================
figure(2)
%contour(xgrid,ygrid,abs(Ezgrid),50);
contourf(xgrid,ygrid,abs(Ezgrid));
xlabel('x (wavelengths)');
ylabel('y (wavelengths)');
axis([-rho rho -rho rho]);
axis square;
colorbar;
%
% Evaluate the exact analytical solution at a boundary half-way between the
% scatterer and the ABC boundary
% =========================================================================
Np=50;
d2p=pi/180;
dist=robj+(rho-robj)/2;
Ezexct=zeros(1,1441); %Ezexct=zeros(1,721);
for I=1:1441 %721
    phi(I)=(I-1)*0.25; %0.5;
    xeval=dist*cos(phi(I)*d2p);
    yeval=dist*sin(phi(I)*d2p);
    for e=1:Nelms     
        x2p=x(elmnod(e,2))-xeval;
        x3p=x(elmnod(e,3))-xeval;
        y2p=y(elmnod(e,2))-yeval;
        y3p=y(elmnod(e,3))-yeval;
        A1=0.5*abs(x2p*y3p-x3p*y2p);
                
        x2p=x(elmnod(e,2))-xeval;
        x1p=x(elmnod(e,1))-xeval;
        y2p=y(elmnod(e,2))-yeval;
        y1p=y(elmnod(e,1))-yeval;
        A2=0.5*abs(x2p*y1p-x1p*y2p);
                
        x1p=x(elmnod(e,1))-xeval;
        x3p=x(elmnod(e,3))-xeval;
        y1p=y(elmnod(e,1))-yeval;
        y3p=y(elmnod(e,3))-yeval;
        A3=0.5*abs(x1p*y3p-x3p*y1p);
                
        x21=x(elmnod(e,2))-x(elmnod(e,1));
        x31=x(elmnod(e,3))-x(elmnod(e,1));
        y21=y(elmnod(e,2))-y(elmnod(e,1));
        y31=y(elmnod(e,3))-y(elmnod(e,1));
        Ae=0.5*(x21*y31-x31*y21);
                
        if abs(Ae-(A1+A2+A3)) < 0.00001*Ae   
            ksi=(y31*(xeval-x(elmnod(e,1)))-x31*(yeval-y(elmnod(e,1))))/(2*Ae);
            ita=(-y21*(xeval-x(elmnod(e,1)))+x21*(yeval-y(elmnod(e,1))))/(2*Ae);
            N1=1-ksi-ita;
            N2=ksi;
            N3=ita;
            Ezeval(I)=N1*Ez(elmnod(e,1))+N2*Ez(elmnod(e,2))+N3*Ez(elmnod(e,3));
        end
    end
    for n=-Np:Np
        Ezexct(I)=Ezexct(I)+E0*j^(-n)*(besselj(n,k0*dist)-besselj(n,k0*robj)/besselh(n,2,k0*robj)*besselh(n,2,k0*dist))*exp(j*n*phi(I)*d2p);
    end
end
%
% Plot the exact analytical solution and the finite element solution along
% a boundary half-way between the scatterer and the ABC boundary
% ========================================================================
figure(3)
plot(phi,abs(Ezexct),'k-',phi,abs(Ezeval),'k--'),legend('Exact','FEM');
xlabel('Observation Angle (degrees)');
ylabel('Electric Field (V/m)');
axis([0 360 0 2*E0]);