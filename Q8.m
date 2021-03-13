
%##########################################################################
%                           Ömer KARAGÖZ 2017                             %
%                         okaragoz@ogu.edu.tr                             %
%                                                                         %
%                   FEM SOLVER USING Q8 MESH ELEMENTS                     %
%This code was written on 01.05.2017 within the scope of METU CE583 course%
%Fully automatically solves for cantilever beam problem                   %
%INSTRUCTOR: Dr. Ozgur KURC                                               %
%MATLAB R2015b - 64bit is used                                            %
%#########################################################################%
% Let's first clear all the previously defined variables
clc; clear all;
% tic and toc used to measure the analyses time
tic;
%-------------------------------------------------------------------------%
%                             User Input                                  %
%-------------------------------------------------------------------------%
%                                                         | endLoad
%  #|                                                     v
%  #|-----------------------------------------------------|
%  #|/////////////////////////////////////////////////////| d
%  #|-----------------------------------------------------|
%  #| <-----                  L                     ------>
E = 25000000;         % Elasticity modulus of the material used
Nu = 0.0;             % Poisson Ratio
d = 1;                % Height of the cantilever beam (m)
t = 1;                % Depth of the cantilever beam (m)
L = 4;                % Length of the cantilever beam (m)
endLoad = -20;        % Load (kN) applied on the free end of the cantilever 
                      % beam endLoad is divided to the nodes on the free
                      % end automatically. endLoad is the total load.
                      
numMeshx = 4;                   % Number of mesh along X axis
numMeshy = 2;                   % Number of mesh along Y axis
%-------------------------------------------------------------------------%
%           !!! Do Not Change Anything Below This Point !!!               %
%-------------------------------------------------------------------------%
xSpacing = (L/numMeshx)/2;     % Length of the mesh element on X axis
ySpacing = (d/numMeshy)/2;     % Length of the mesh element on Y axis
numElem = numMeshx*numMeshy;       % Number of total elements
% Number of nodes
numNode = (numMeshx*2+1)*(numMeshy+1)+numMeshy*(numMeshx+1);
%-------------------------------------------------------------------------%
%                        Coordinate Calculation                           %
%-------------------------------------------------------------------------%
% This part calculates the coordinates for all the nodes (1 through to end)
nodeCoord = zeros(numNode,2);
c = 0;
for i=1:(numMeshy*2+1)
    if mod(i,2) == 1
        for j=1:(2*numMeshx+1)
            c = c+1;
            xCoord = xSpacing*(j-1);
            yCoord = ySpacing*(i-1);
            nodeCoord(c,:) = [xCoord yCoord];
        end
    else
        for j=1:(numMeshx+1)
            c = c+1;
            xCoord = (2*xSpacing)*(j-1);
            yCoord = ySpacing*(i-1);
            nodeCoord(c,:) = [xCoord yCoord];
        end
    end
end
%-------------------------------------------------------------------------%
%                    Connectivity Matrix Computation                      %
%-------------------------------------------------------------------------%
% This part calculates the connectivity matrix. It can calculate any number
% of mesh count along X and Y directions automatically. 
conn = zeros(numElem,8);
rowdiff = numMeshx*2+numMeshx+2;
rowdiffsmall = numMeshx*2+1;
c = 0;
for j=0:numMeshy-1
    for i=0:numMeshx-1
        c=c+1;
        c1 = (2*i+1)+rowdiff*j;
        c2 = c1+2;
        c4 = c1+rowdiff;
        c3 = c4+2;
        c5 = c1+1; 
        c7 = c4+1;
        c8 = c1+rowdiffsmall-i;
        c6 = c8+1;
        conn(c,:) = [c1; c2; c3; c4; c5; c6; c7; c8];
    end
end
%-------------------------------------------------------------------------%
%                        Load Matrix Computation                          %
%-------------------------------------------------------------------------%
% This part divides the endLoad to the nodes which are on the free end.
load = zeros(numMeshy*2+1,3);
load(1,:) = [rowdiffsmall 0 endLoad/(numMeshy*2+1)];
for i=1:numMeshy*2
    load(i+1,:) = [load(i,1)+round((numMeshx+1)*(abs(cos((i+1)*pi()/2)))...
           +rowdiffsmall*(abs(cos(i*pi()/2))),0) 0 endLoad/(numMeshy*2+1)];
end
numLoad = size(load,1);         % Number of loads
%-------------------------------------------------------------------------%
%                      Support Matrix Computation                         %
%-------------------------------------------------------------------------%
% This part calculates the support matrix. First it finds the nodes which
% are on the fixed end and than fixes that nodes.
s = zeros(numMeshy*2+1,3);
s(1,:) = [1 1 1];
for i=1:numMeshy*2                    % s=[nodeID Ux Uy] 1-> Fixed 0-> Free
    s(i+1,:) = [load(i,1)+1 1 1];
end
numSupport = size(s,1);               % Number of supports
c=0;
for i=1:numElem
    nodeCoordElem = zeros(8,2);
    for j=1:8
        nodeCoordElem(j,:) = nodeCoord(conn(i,j),:);
    end
    c=c+1;
    nodeCoordElemList(c).nodeCoordElem = nodeCoordElem;
end
EM = E/(1-Nu^2).*[1 Nu 0; Nu 1 0; 0 0 (1-Nu)/2];      % Constitutive Matrix
syms ksi ita
% Serendipity Shape Functions
n1 = 0.25*(-(-1+ksi)*(-1+ita)*(ksi+ita+1));
n2 = 0.25*(1+ksi)*(-1+ita)*(-ksi+ita+1);
n3 = 0.25*(1+ksi)*(1+ita)*(ksi+ita-1);
n4 = 0.25*(-(-1+ksi)*(1+ita)*(-ksi+ita-1));
n5 = 0.5*(-1+ksi^2)*(-1+ita);
n6 = -0.5*(1+ksi)*(-1+ita^2);
n7 = -0.5*(-1+ksi^2)*(1+ita);
n8 = 0.5*(-1+ksi)*(-1+ita^2);
N8 = [n1; n2; n3; n4; n5; n6; n7; n8];
% Differantion of the shape functions calculated to use later.
dnksi = [diff(n1,ksi) diff(n2,ksi) diff(n3,ksi) diff(n4,ksi)...
         diff(n5,ksi) diff(n6,ksi) diff(n7,ksi) diff(n8,ksi)];
dnita = [diff(n1,ita) diff(n2,ita) diff(n3,ita) diff(n4,ita)...
         diff(n5,ita) diff(n6,ita) diff(n7,ita) diff(n8,ita)];
D8 =[dnksi
     dnita];
%-------------------------------------------------------------------------%
%                 Element Stiffness Matrix Calculation                    %
%-------------------------------------------------------------------------%
% This part calculates the element stiffness matrices for each element.
parfor g=1:numElem
    jac=D8*nodeCoordElemList(g).nodeCoordElem;    % Jacobian Matrix
    capGamma = inv(jac);                          % Inverse Jacobian Matrix
    B = sym(zeros(3,16));   % B: Strain-Displacement relation matrix (3x16)
    for k=0:7
        m = 2*k+1;
        n = m+1;
        B(1,m) = capGamma(1,1)*dnksi(k+1)+capGamma(1,2)*dnita(k+1);
        B(2,n) = capGamma(2,1)*dnksi(k+1)+capGamma(2,2)*dnita(k+1);
        B(3,m) = B(2,n);
        B(3,n) = B(1,m);
    end
    kElemList(g).kElem = int(int(transpose(B)*EM*B*t*det(jac),ksi,-1,1)...
                                                                ,ita,-1,1);
    kElemList(g).B = B;
end
%-------------------------------------------------------------------------%
%                  Global Stiffness Matrix Calculation                    %
%-------------------------------------------------------------------------%
% This part calculates the global stiffness matrix. Basically; for each
% element it takes 4x4 part from the element stiffness matrix and puts to
% the correct spot on the global stiffness matrix. This process loops until
% all elements all parts placed in to the global stiffness matrix.
kSystem = zeros(numNode*2,numNode*2);
for k=1:numElem
    kElem = kElemList(k).kElem;
    for j=1:8
        for i=1:8
            n = conn(k,i);
            m = conn(k,j);
            kSystem(2*n-1,2*m-1) = kSystem(2*n-1,2*m-1)+kElem(2*i-1,2*j-1);
            kSystem(2*n-1,2*m) = kSystem(2*n-1,2*m)+kElem(2*i-1,2*j);
            kSystem(2*n,2*m-1) = kSystem(2*n,2*m-1)+kElem(2*i,2*j-1);
            kSystem(2*n,2*m) = kSystem(2*n,2*m)+kElem(2*i,2*j);
        end
    end
end
% kSystem is copied before the support conditions were applied for later
% use
kSystemWOsup = kSystem;
%-------------------------------------------------------------------------%
%           Apply Support Conditions to the Global Stiffnes Matrix        %
%-------------------------------------------------------------------------%
% This part makes zeros all columns and rows where the supports are except
% the diagonal element of the matrix. Diagonal element set to 1. I choose
% this method because with this way sort of displacement evaulated are not
% changes. And later it is easier to use evaluated values. Only negative
% side of this approach is we have to be careful not to put force where the
% support is fixed.
for i=1:numSupport
    n = s(i,1);
    if (s(i,2) == 1)
        kSystem(2*n-1,:) = 0;
        kSystem(:,2*n-1) = 0;
        kSystem(2*n-1,2*n-1) = 1;
    end
    if (s(i,3) == 1)
        kSystem(2*n,:) = 0;
        kSystem(:,2*n) = 0;
        kSystem(2*n,2*n) = 1;
    end
end
%-------------------------------------------------------------------------%
%                       Load Vector Computation                           %
%-------------------------------------------------------------------------%
% In this part load vector created. If there is a load vector get the value
% from load matrix. If not not load value set to zero.
loadVector = zeros(numNode*2,1);
for i=1:numLoad
    n = load(i,1);
    loadVector(2*n-1) = load(i,2);
    loadVector(2*n) = load(i,3);
end
%-------------------------------------------------------------------------%
%                       Displacement Computation                          %
%-------------------------------------------------------------------------%
u = kSystem\loadVector;
for i=1:numNode                             % Displacements for each node
    uNodeList(i).Node = [u(2*i-1); u(2*i)];
end
for i=1:numElem                             % Displacements for each elemet
    temp = [];
    for j=1:8
        temp = [temp; uNodeList(conn(i,j)).Node];
    end
    uElemList(i).Elem = temp;
end
disp('# Tip Displacements\n');
for i=1:(numMeshy*2+1)          % Display all displacements on the free end
    fprintf('%10.7f \n',u(load(i,1)*2));
end
% Displacement from beam theory (From Vertical Load + From Shear)
BeamTheory = (endLoad*L^3)/(3*E*(t*d^3/12))+...
                                      (endLoad*L)/((E/2*(1-Nu))*(5/6)*d*t);
fprintf('Beam Theory: %10.7f \n',BeamTheory);
%-------------------------------------------------------------------------%
%                         Stress Computation                              %
%-------------------------------------------------------------------------%
SCP = [-1 -1,; 1 -1; 1 1; -1 1];                % Stress calculation points
   
parfor i=1:numElem
   uelem = uElemList(i).Elem';
   temp = [];
   fElemList(i).fElem = double(kElemList(i).kElem*uelem');
   sigma = EM*kElemList(i).B*uelem';
   temp = [];
   for j=1:4
       temp = [temp double(subs(subs(sigma,ksi,SCP(j,1)),ita,SCP(j,2)))]; 
   end
   sigmaElemList(i).sigma = temp;
end
%-------------------------------------------------------------------------%
%                          Force Computation                              %
%-------------------------------------------------------------------------%
fSystem = kSystemWOsup*u;
% End of the analysis
toc;
drawSystem;                % Optional step. This scripts draws deformed and
                           % undeformed shape with colored stress values.
