function eta=SteeringWNR(SIGax)

dS=size(SIGax);
dim=dS(1);
Oa=dS(3);
Ia=dS(4);
Daxl=Dax_MATRIX(Oa,Ia);
nL=Oa^Ia;  %Count the deterministic vertices

cvx_begin SDP 
variable eta
expression SIGax_Noise(dim,dim,Oa,Ia)
variable sig_l(dim,dim,nL) semidefinite complex
expression sumE_l(dim,dim,Oa,Ia)
for a=1:Oa
    for x=1:Ia
        SIGax_Noise(:,:,a,x)=eta*SIGax(:,:,a,x)+(1-eta)*eye(dim)/Oa/dim;
        for l=1:nL
            sumE_l(:,:,a,x) = sumE_l(:,:,a,x) + Daxl(a,x,l)*sig_l(:,:,l);
        end
        SIGax_Noise(:,:,a,x)==sumE_l(:,:,a,x);
    end
end

maximise eta

cvx_end

JM=eta;
end
