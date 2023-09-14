%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Requires: QETLAB, cvx
%Last update: 03/05/2023

%This script considers that the state experimentally obtained via which is close to a qubit isotropic state with visibility alpha_target=0.6373;
%The goal of this script is to prove that the experimental state admits an LHV model for arbitrary projective measurements

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

cvx_begin SDP
variable rho_PPT(d^2 d^2) complex semidefinite %Arbitrary separable two-qubit state (unnormalised)
variable map_CP1(d^2 d^2) complex semidefinite %Choi operator of a CP map (unnormalised)
variable map_CP2(d^2 d^2) complex semidefinite %Choi operator of a CP map (unnormalised)
variable eta %critical visibility
expression map %Choi operator of a positive map (unnormalised)
expression rho_local %IsotripicLocal after an arbitrary positive map in its first component

PartialTranspose(rho_PPT)>=0; %rho_PPT is PPT
map=map_CP1+PartialTranspose(map_CP2); %Canonical form of the choi operaotor of an arbitrary positive
PartialTrace(map,2,[d d])==trace(map)*eye(d)/d; %Ensures that map is proportional to a trace-preserving map
rho_local=PartialTrace(kron(map,eye(d))*kron(eye(d),PartialTranspose(IsotropicLocal,2,[d d])),2,[d d d]); %Applying linear maps via Choi isomorphism. rho_local is an unnormaised state which necessarily has a local model

eta*rho_target+(1-eta)*eye(d^2)/d^2 == rho_local + rho_PPT; %Write rho_target as a convex combination of states with LHV model
% eta<=1; %Condition to prevent eta to be greater than 1 (not really needed, but convenient, just to simplify analysis)

maximise eta;
cvx_end

eta=eta

if eta>=1-10^(-8) %Leave some error margin for the SDP
	display('The target state could be written as a convex combination of states with LHV, hence it has a LHV model.')
end

% The code below simply tests the fidelity betwen rho_target and the isotropic state with visibility alpha_target=0.6373 
% alpha_target=0.6373;
% display(' ')
% display('The fidelity between the target state and an isotropic state with visibility alpha_target=0.6373 is: ')
% experimetal_fidelity=Fidelity(rho_target,alpha_target*PhiP+(1-alpha_target)*eye(d^2)/d^2)

%%%%%%%% This part below is not needed, but it might be used to double-check that everything behaves as it should (unncomment lines to obtain the results)
%For the sake of concreteness, below we explicitly decompose rho_target as a convex comination of a state with LHV and a PPT state. 
rho_PPT_normalised=rho_PPT/trace(rho_PPT);
q=trace(rho_local);
% eig(map_CP1)
% eig(map_CP2)
% eig(rho_PPT)
% eig(PartialTranspose(rho_PPT))
map_normalised=map/trace(map)*d;
rho_local_normalised=PartialTrace(kron(map_normalised,eye(d))*kron(eye(d),PartialTranspose(IsotropicLocal,2,[d d])),2,[d d d]); %Applying linear maps via Choi isomorphism. rho_local is an unnormaised state which necessarily has a local model
rho_LHV=q*rho_local_normalised+(1-q)*rho_PPT_normalised;
%norm(rho_target-rho_LHV)