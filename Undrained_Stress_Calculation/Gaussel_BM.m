% this function solve the banded linear system using Gaussian elimination

function U = Gaussel_BM (K_bm,Residual)

U = zeros(length(Residual),1);
[row,col] = size(K_bm);
ma = 0.5*(col+1) - 1;


%-------------------------%
%   FORWARD ELIMINATION   %
%-------------------------%

NN1 = row - 1;

for IQ=1:NN1
    
    I1 = IQ + 1;
    
    for i=I1:row
        if (ma+1+IQ-i) > 0
            temp = K_bm(i,ma+1+IQ-i)/K_bm(IQ,ma+1);
            Residual(i,1)=Residual(i,1)-Residual(IQ,1)*temp;
            for j=col:-1:ma+1
                j1=j+IQ-i;
                if j1>0
                    K_bm(i,j1) = K_bm(i,j1) - K_bm(IQ,j)*temp;
                end
            end
        end
    end
end


%---------------------%
%   BACK SUBSTITUTE   %
%---------------------%

if K_bm(row,ma+1)~=0
    U(row) = Residual(row,1)/K_bm(row,ma+1);
else
    U(row) = 0;
end

I1 = ma + 1;

for I2=1:NN1
    IQ = row - I2;
    if IQ<(row-ma)
        for i=1:ma
            Residual(IQ,1) =  Residual(IQ,1)-K_bm(IQ,I1+i)*U(IQ+i,1);
        end
    else
        for i=1:I2
            Residual(IQ,1) =  Residual(IQ,1)-K_bm(IQ,I1+i)*U(IQ+i,1);
        end
    end
    U(IQ,1)=Residual(IQ,1)/K_bm(IQ,ma+1);
end
end