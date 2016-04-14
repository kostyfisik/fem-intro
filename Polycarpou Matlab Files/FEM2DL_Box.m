%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NAME OF FILE: FEM2DL_Box.m
%
% PURPOSE: This Matlab code solves a two-dimensional boundary-value problem
% using the finite element method. The boundary-value problem at hand
% corresponds to the Laplace's equation as applied to a rectangular domain 
% characterized by a simple medium (isotropic, linear, and homogeneous)
% with dielectric constant epsr. The finite element method is based on sub-
% dividing the two-dimensional domain into Ne linear triangular elements. 
% Dirichlet boundary conditions are applied to all four metallic walls: 
%
% V=0 on the left sidewall 
% V=0 on the right sidewall
% V=0 on the bottom sidewall
% V=V0 on the top wall which is separated by tiny gaps from the two sidewalls
%
% The dimensions of the domain are WxH, where W=Width and H=Hight. 
% The user has the freedom to modify any of the defined input variables 
% including domain width and height (W & L), top wall voltage (V0), and number 
% of bricks along the x- and y-axes. Note that each brick is then subdivided 
% into two triangles.
%
% All quantities are expressed in the SI system of units.
%
% The Primary Unknown Quantity is the Electric Potential.
% The Secondary Unknown Quantity is the Electric Field.
%
% Written by Anastasis Polycarpou (Last updated: 8/12/2005)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear % Clear all variables

% Define Rectangular Geometry
% ===========================
W=1;
H=1;

% Define the number of bricks in the x and y directions
% =====================================================
XBRICKS=20;
YBRICKS=20;

% Define the Voltage (Electric potential) on the top wall
% =======================================================
V0=1;

% Generate the triangular mesh
% ============================
XLLC=0;  % x-coordinate of the Lower Left Corner
YLLC=0;  % y-coordinate of the Lower Left Corner
XURC=W;  % x-coordinate of the Upper Right Corner
YURC=H;  % y-coordinate of the Upper Right Corner
TBRICKS=XBRICKS*YBRICKS;   % Total number of bricks
Dx=(XURC-XLLC)/XBRICKS;  % Edge length along x-direction
Dy=(YURC-YLLC)/YBRICKS;  % Edge length along y-direction

TNNDS=(XBRICKS+1)*(YBRICKS+1);    % Total number of nodes
TNELMS=2*TBRICKS;    % Total number of triangular elements
%
% Print the number of nodes and the number of elements
% ====================================================
fprintf('The number of nodes is: %6i\n',TNNDS)
fprintf('The number of elements is: %6i\n',TNELMS)

for I=1:TNNDS
    X(I)=mod(I-1,XBRICKS+1)*Dx;
    Y(I)=floor((I-1)/(XBRICKS+1))*Dy;
end

% Connectivity Information
% ========================
for I=1:TBRICKS
    ELMNOD(2*I-1,1)=I+floor((I-1)/XBRICKS);
    ELMNOD(2*I-1,2)=ELMNOD(2*I-1,1)+1+(XBRICKS+1);
    ELMNOD(2*I-1,3)=ELMNOD(2*I-1,1)+1+(XBRICKS+1)-1;
    
    ELMNOD(2*I,1)=I+floor((I-1)/XBRICKS);
    ELMNOD(2*I,2)=ELMNOD(2*I,1)+1;
    ELMNOD(2*I,3)=ELMNOD(2*I,1)+1+(XBRICKS+1);
end
%
% Plot the mesh
% =============
figure(1)
triplot(ELMNOD,X,Y);

% Define the constants for each element
% =====================================
for e=1:TNELMS
    alphax(e)=1;
    alphay(e)=1;
    beta(e)=0;
    g(e)=0;
end

% Definition of Dirichlet Boundary Conditions
% ===========================================
TNEBC=0;
for I=1:TNNDS
    if X(I) == XLLC || X(I) == XURC || Y(I) == YLLC
        TNEBC=TNEBC+1;
        EBCNOD(TNEBC)=I;
        EBCVAL(TNEBC)=0;
    elseif Y(I) == YURC
        TNEBC=TNEBC+1;
        EBCNOD(TNEBC)=I;
        EBCVAL(TNEBC)=V0;
    end
end

% Definition of Mixed Boundary Conditions
% =======================================
TNMBC=0;

% Initialization of the global K matrix and right-hand side vector
% ================================================================
K=sparse(TNNDS,TNNDS);
b=zeros(TNNDS,1);

% Form the element matrices and assemble to the global matrix
% ===========================================================
for e=1:TNELMS
    x21=X(ELMNOD(e,2))-X(ELMNOD(e,1));
    x31=X(ELMNOD(e,3))-X(ELMNOD(e,1));
    x32=X(ELMNOD(e,3))-X(ELMNOD(e,2));
    x13=X(ELMNOD(e,1))-X(ELMNOD(e,3));
    y12=Y(ELMNOD(e,1))-Y(ELMNOD(e,2));
    y21=Y(ELMNOD(e,2))-Y(ELMNOD(e,1));
    y31=Y(ELMNOD(e,3))-Y(ELMNOD(e,1));
    y23=Y(ELMNOD(e,2))-Y(ELMNOD(e,3));
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
        % There is no boundary Gamma2 !!!!
    
    % Assemble element matrices & vectors into the global K matrix and b
    % vector
    % ==================================================================
    for I=1:3
        for J=1:3
            K(ELMNOD(e,I),ELMNOD(e,J))=K(ELMNOD(e,I),ELMNOD(e,J))+Ke(I,J);
        end
        b(ELMNOD(e,I))=b(ELMNOD(e,I))+ge(I);
    end
end

% Imposition of Dirichlet boundary conditions
% ===========================================
for I=1:TNEBC
    for J=1:TNNDS
        if J ~= EBCNOD(I) 
            b(J)=b(J)-K(J,EBCNOD(I))*EBCVAL(I);
        end
    end
    K(:,EBCNOD(I))=0;
    K(EBCNOD(I),:)=0;
    K(EBCNOD(I),EBCNOD(I))=1;
    b(EBCNOD(I))=EBCVAL(I);
end

% Solution of the global matrix system
% ====================================
V=K\b;

% Generate the solution over a grid and plot it
% =============================================
[Xgrid,Ygrid]=meshgrid(XLLC:0.01*(XURC-XLLC):XURC,YLLC:0.01*(YURC-YLLC):YURC);
Vgrid=zeros(101,101);
for I=1:101
    for J=1:101
        for e=1:TNELMS
         
            x2p=X(ELMNOD(e,2))-Xgrid(I,J);
            x3p=X(ELMNOD(e,3))-Xgrid(I,J);
            y2p=Y(ELMNOD(e,2))-Ygrid(I,J);
            y3p=Y(ELMNOD(e,3))-Ygrid(I,J);
            A1=0.5*abs(x2p*y3p-x3p*y2p);
                
            x2p=X(ELMNOD(e,2))-Xgrid(I,J);
            x1p=X(ELMNOD(e,1))-Xgrid(I,J);
            y2p=Y(ELMNOD(e,2))-Ygrid(I,J);
            y1p=Y(ELMNOD(e,1))-Ygrid(I,J);
            A2=0.5*abs(x2p*y1p-x1p*y2p);
                
            x1p=X(ELMNOD(e,1))-Xgrid(I,J);
            x3p=X(ELMNOD(e,3))-Xgrid(I,J);
            y1p=Y(ELMNOD(e,1))-Ygrid(I,J);
            y3p=Y(ELMNOD(e,3))-Ygrid(I,J);
            A3=0.5*abs(x1p*y3p-x3p*y1p);
                
            x21=X(ELMNOD(e,2))-X(ELMNOD(e,1));
            x31=X(ELMNOD(e,3))-X(ELMNOD(e,1));
            y21=Y(ELMNOD(e,2))-Y(ELMNOD(e,1));
            y31=Y(ELMNOD(e,3))-Y(ELMNOD(e,1));
            Ae=0.5*(x21*y31-x31*y21);
                
            if abs(Ae-(A1+A2+A3)) < 0.00001*Ae   
                 ksi=(y31*(Xgrid(I,J)-X(ELMNOD(e,1)))-x31*(Ygrid(I,J)-Y(ELMNOD(e,1))))/(2*Ae);
                 ita=(-y21*(Xgrid(I,J)-X(ELMNOD(e,1)))+x21*(Ygrid(I,J)-Y(ELMNOD(e,1))))/(2*Ae);
                 N1=1-ksi-ita;
                 N2=ksi;
                 N3=ita;
                 Vgrid(I,J)=N1*V(ELMNOD(e,1))+N2*V(ELMNOD(e,2))+N3*V(ELMNOD(e,3));
            end
        end
    end
end
 
% Plot the finite element solution of V using a contour plot
% ==========================================================
figure(2)
contourf(Xgrid,Ygrid,Vgrid,15);
xlabel('x');
ylabel('y');
colorbar;

% Exact Analytical Solution
% =========================
N=100;
Vexact=zeros(101,101);
for i=1:101
    for j=1:101
        for k=1:N
            Update=(4/pi)*sin((2*k-1)*pi/W*Xgrid(i,j))*sinh((2*k-1)*pi/W*Ygrid(i,j))/((2*k-1)*sinh((2*k-1)*pi*H/W));
            if abs(Update) > abs(0.001*Vexact(i,j))
                Vexact(i,j)=Vexact(i,j)+Update;
            end
        end
    end
end

% Error based on L2 Norm
% ======================
L2error=0;
for i=1:101
    for j=1:101
        L2error=L2error+(Vexact(i,j)-Vgrid(i,j))^2;
    end
end
L2error=sqrt(L2error/(101*101))

    
