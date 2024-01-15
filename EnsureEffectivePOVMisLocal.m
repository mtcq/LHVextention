%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Requires: QETLAB, cvx
%Last update: 15/01/2024

%The goal of this script is to ensure that the effective 4-outcome, 4-input
%POVM perfomed by B and C after the isometric channel cannot lead to Bell
%nonlocality. For this goal, we show that this set of effective POVMs with
%the state from our experiment is unsteerable. We also show that the white
%noisy robustness of this effective POVM is 0.7746, hence considerably
%grater than p_Designolle = 0.6875, the visibility for which the two-qubit
%Werner state is ensure to be Bell local.

clear all
d=2; %Dimension of each subsystem
%We define the tomographed state by writing is real and imaginary parts
rho_real=[0.39433196722500463    0.004063527194669355    -0.018455908278113917    0.31320663947255156;
    0.004063527194669355    0.08979326112974952    -0.0006089441954363349    0.01382540917454235;
    -0.018455908278113917    -0.0006089441954363349    0.08886298939962693    -0.0050651520542264504;
    0.31320663947255156    0.01382540917454235    -0.0050651520542264504    0.42701178224561886];

rho_imag=[0.0    0.011814846261924255    -0.0035099174111995275    -0.034383852546323104;
    -0.011814846261924255    0.0    -0.010418163316369543    -0.012206421234440292;
    0.0035099174111995275    0.010418163316369543    0.0    -0.012240162478548458;
    0.034383852546323104    0.012206421234440292    0.012240162478548458    0.0];

rho_target=rho_real+1i*rho_imag;
%We now define the best known isotropic state has an LHV model
ketPhiP=[1; 0; 0; 1;]/sqrt(2); %Defines the maximally entangled state
PhiP=ketPhiP*ketPhiP'; %Define its density matrix
p_Designolle=0.6875; %Visibility for which the two-qubit isotropic state is known to have an LHV model for projective measurements (https://arxiv.org/abs/2302.04721)
IsotropicLocal=p_Designolle*PhiP+(1-p_Designolle)*eye(4)/4; %Best known isotropic state has an LHV model

X=[0 1;1 0];
Y=[0 -sqrt(-1);sqrt(-1) 0];
Z=[1 0;0 -1];

%Define Observables for A, B, and C
A0=(-X-Z)/sqrt(2);
A1=(X-Z)/sqrt(2);
A2=-Y;

B0=(sqrt(2)*X+Y)/sqrt(3);
B1=(sqrt(2)*X-Y)/sqrt(3);
C0=Z;
C1=X;

%Define POVM elements for A, B, and C
A(:,:,1,1)=(A0+eye(2))/2;
A(:,:,1,2)=(A1+eye(2))/2;
A(:,:,1,3)=(A2+eye(2))/2;
A(:,:,2,1)=eye(2)-A(:,:,1,1);
A(:,:,2,2)=eye(2)-A(:,:,1,2);
A(:,:,2,3)=eye(2)-A(:,:,1,3);

B(:,:,1,1)=(B0+eye(2))/2;
B(:,:,1,2)=(B1+eye(2))/2;
B(:,:,2,1)=eye(2)-B(:,:,1,1);
B(:,:,2,2)=eye(2)-B(:,:,1,2);

C(:,:,1,1)=(C0+eye(2))/2;
C(:,:,1,2)=(C1+eye(2))/2;
C(:,:,2,1)=eye(2)-C(:,:,1,1);
C(:,:,2,2)=eye(2)-C(:,:,1,2);

%Define the Isometry V= 1/sqrt(2) ( (|00>-|11>)<0| +  (-|01>-|10>)<1| ) 
V=[1;0;0;-1]*[1 0]/sqrt(2)+[0;-1;-1;0]*[0 1]/sqrt(2);

%Define the tensor produc of POVMs B and C
BC(:,:,1,1)=kron(B(:,:,1,1),C(:,:,1,1));
BC(:,:,2,1)=kron(B(:,:,1,1),C(:,:,2,1));
BC(:,:,3,1)=kron(B(:,:,2,1),C(:,:,1,1));
BC(:,:,4,1)=kron(B(:,:,2,1),C(:,:,2,1));

BC(:,:,1,2)=kron(B(:,:,1,1),C(:,:,1,2));
BC(:,:,2,2)=kron(B(:,:,1,1),C(:,:,2,2));
BC(:,:,3,2)=kron(B(:,:,2,1),C(:,:,1,2));
BC(:,:,4,2)=kron(B(:,:,2,1),C(:,:,2,2));

BC(:,:,1,3)=kron(B(:,:,1,2),C(:,:,1,1));
BC(:,:,2,3)=kron(B(:,:,1,2),C(:,:,2,1));
BC(:,:,3,3)=kron(B(:,:,2,2),C(:,:,1,1));
BC(:,:,4,3)=kron(B(:,:,2,2),C(:,:,2,1));

BC(:,:,1,4)=kron(B(:,:,1,2),C(:,:,1,2));
BC(:,:,2,4)=kron(B(:,:,1,2),C(:,:,2,2));
BC(:,:,3,4)=kron(B(:,:,2,2),C(:,:,1,2));
BC(:,:,4,4)=kron(B(:,:,2,2),C(:,:,2,2));

%Defines the effective POVM BC after the isometry V
for a=1:4
    for x=1:4
        BCV(:,:,a,x)=V'*BC(:,:,a,x)*V;
    end
end

for x=1:4
    for a=1:4
        SIGax(:,:,a,x)=PartialTrace(rho_target*kron(BCV(:,:,a,x),eye(2)),1);
    end
end

eta_critical = SteeringWNR(SIGax);
WNR_JM=JM(BCV);
display('The robustness to white noise of the assemblage with rho_target is: ')
eta_critical=eta_critical

display('The robustness to white noise of the 4-outcome effective POVMs is: ')
WNR_JM=WNR_JM

