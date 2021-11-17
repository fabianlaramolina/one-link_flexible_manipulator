%Paper Title: Stochastic Finite-Element Modeling of Flexible Manipulators 
%Dynamic Model of the One-Link Flexible Manipulator by using Finite Element
%Method
%Fabian Andres Lara-Molina
%Federal University of Technology - Paraná
%Lara.f8@gmail.com
%APPENDIX A
%2021

function One_Link_Flexible_Manipulator_1_DOF_Appendix_A()
%------------------------------------------------------
%------------------------------------------------------
%Definition of the variables:
% rho  - mass density of the link
% A - cross-section area of the link
% E - Young's Modulus of the link
% I - Moment of cross section area - link
% l - link length
% y - r=[x y]' y coordinate of the element Cartesian position 
% dy - time derivative of y
% u3 - Translational elastic degee of freedom
% u4 - Rotational elastic degee of freedom 
% du3 - time derivative of u3
% du4 - time derivative of u4
% phi_1 - Joint angle
% dphi_1 - time derivative of phi_1

%psi_1=[u3 u4] - Elastic degrees of freedom
%q_1=[phi_1 u3 u4] - Manipulator Degrees of freedom
%dq_1=[dphi_1 du3 du4] - time derivative of q_1
%------------------------------------------------------
%------------------------------------------------------

clear all
clc
syms rho A E I
syms l y dy 
syms u3 u4 du3 du4 phi_1 dphi_1 

j=1;                    %Number of elements
ne=4;                   %Number of nodes of the element
nd=(j+1)*ne/2+1;        %Number of nodes of the beam

%==========================================================================
%1), 2) and 3)Matrices for the Total Energy

%Elementary Matrices Definition----------------------------------------------------------
K=Stiffness_Matrix(l,E,I); %Stiffness Massa Matriz Definition
M=Mass_Matrix(l,y,rho,A);  %Elementary Massa Matriz for j=1
psi_1=[0 0 u3 u4]; %Elastic degees of freedom of the element
                   %Constraints:u1=0 u2=0 - Link is atached to the hub

M=Sub_non_linear_term(psi_1,M,y,l); %Function to include non-linear terms
                                    %within the mass elementary matrix

%==========================================================================
% 4) Boundary conditions:
%Constraints - flexural and rotational displacements at the hub as zero
%Stiffness Matrix: Constraints u1=0 u2=0
K_T_nom=sym(zeros(nd-2,nd-2));
K_T_nom(2:nd-2,2:nd-2)=K(4:nd,4:nd);
K_T_nom(2:nd-2,1)=K(4:nd,1);
K_T_nom(1,:)=[K(1,1) K(1,4:nd)];

%Mass Matrix: Constraints u1=0 u2=0
M_T_nom=sym(zeros(nd-2,nd-2));
M_T_nom(2:nd-2,2:nd-2)=M(4:nd,4:nd);
M_T_nom(2:nd-2,1)=M(4:nd,1);
M_T_nom(1,:)=[M(1,1) M(1,4:nd)];



%==========================================================================
% 5) Lagrange Formulation
%Coriolis/centripetal vector
%C=dM(q)*dq-1/2*d/dq*(dq*M(q)*dq) - Eq(8)

dq_1=[dphi_1 du3 du4];

M1=simplify(dq_1*M_T_nom*dq_1.');

h1=-1/2*[diff(M1,phi_1);diff(M1,u3);diff(M1,u4)];

dM=fulldiff(M_T_nom,{phi_1,u3,u4});

h2=simplify(dM*[dphi_1;du3;du4]);
%Outputs: 
M_T_nom  %Total inertia Matrix
h=simplify(h1+h2) %Coriolis/Centripetal vector
K_T_nom %Total Stiffness Matrix


function M_i=Mass_Matrix(l,y,rho,A)
M_i=(rho*A*l)/420*[120*l^2+y^2 63*l 14*l^2 147*l -21*l^2;...
    63*l 156 22*l 54 -13*l;...
    14*l^2 22*l 4*l^2 13*l -3*l;...
    147*l 54 13*l 156 -22*l;...
    -21*l^2 -13*l -3*l^2 -22*l 4*l^2];


function K_i=Stiffness_Matrix(l,E,I)
K_i=(E*I)/l^3*[0 0 0 0 0;...
               0 12 6*l -12 6*l;...
               0 6*l 4*l^2 -6*l 2*l^2;...
               0 -12 -6*l 12 -6*l;...
               0 6*l  2*l^2 -6*l 4*l^2];

function M_i=Sub_non_linear_term(psi_1,M_1,y,l)

P=l/420*[156 22*l 54 -13*l;22*l 4*l^2 13*l -3*l^2;54 13*l 156 -22*l;-13*l -3*l^2 -22*l 4*l^2];
yy=simplify(psi_1*P*psi_1.');
M_i=subs(M_1,{y^2},{yy});


