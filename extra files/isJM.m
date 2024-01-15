function isJM=isJM(Max)
%Max(:,:,Oa,Ia)

dS=size(Max);
dim=dS(1);
Oa=dS(3);
Ia=dS(4);
Daxl=Dax_MATRIX(Oa,Ia);
nL=Oa^Ia;  %Count the deterministic vertices

cvx_begin SDP

variable E_l(dim,dim,nL) semidefinite complex
expression sumE_l(dim,dim,Oa,Ia)

for a=1:Oa
    for x=1:Ia
        for l=1:nL
            sumE_l(:,:,a,x) = sumE_l(:,:,a,x) + Daxl(a,x,l)*E_l(:,:,l);
        end
        Max(:,:,a,x)==sumE_l(:,:,a,x);
    end
end
cvx_end

end

