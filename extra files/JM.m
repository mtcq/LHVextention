function JM=JM(Max)
%Function that Evaluetes the white noisy robstness of a set of POVMs Max(:,:,a,x)
dS=size(Max);
dim=dS(1);
Oa=dS(3);
Ia=dS(4);

cvx_begin SDP
    variable eta
    expression Max_Noise(dim,dim,Oa,Ia)
for a=1:Oa
    for x=1:Ia
        Max_Noise(:,:,a,x)=eta*Max(:,:,a,x)+(1-eta)*eye(dim)/Oa;
    end
end
     maximise eta
     isJM(Max_Noise);
   
cvx_end

JM=eta;
end
