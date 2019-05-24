function NMI = community_nmi(A,B)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = numel(A);
Ca = max(A);
Cb = max(B);
N1 = zeros(Ca*Cb,1);
D1 = zeros(Ca,1);
D2 = zeros(Cb,1);

n1_k = 0;

for i=1:Ca
    N_idot=sum(A==i);
    D1(i)=N_idot*log(N_idot/N);
    for j=1:Cb
        n1_k=n1_k+1;
        N_dotj = sum(B==j);
        N_ij=sum(((A==i)+(B==j))==2);
        if N_ij==0
            N1(n1_k)=0;
        else
            N1(n1_k)=N_ij*log((N_ij*N)/(N_idot*N_dotj));
        end
        D2(j)=N_dotj*log(N_dotj/N);  
    end
end

NMI = (-2*sum(N1)) / (sum(D1) + sum(D2));

end

