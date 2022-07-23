% this function solve a linear system using Gaussian elimination

function U = gaussel (Kglobal,Residual)

U = zeros(length(Residual),1);
[row,col] = size(Kglobal);
%foward elimination
NN1 = row - 1;
for IQ=1:NN1
    I1 = IQ + 1;
    for i=I1:row
        Residual(i,1)=Residual(i,1)-Kglobal(i,IQ)*Residual(IQ)/Kglobal(IQ,IQ);
        for j=col:-1:IQ
            Kglobal(i,j)=Kglobal(i,j)-Kglobal(IQ,j)*(Kglobal(i,IQ)/Kglobal(IQ,IQ));
            IQ;
            i;
        end
    end
end

%back substitute
if Kglobal(row,col)~=0
    U(row) = Residual(row,1)/Kglobal(row,col);
else
    U(row) = 0;
end
for I2=1:NN1
    IQ = row - I2;
    I1 = IQ + 1;
    for i=I1:row
        Residual(IQ,1) =  Residual(IQ,1)-Kglobal(IQ,i)*U(i,1);
    end
    U(IQ)=Residual(IQ)/Kglobal(IQ,IQ);
end
end