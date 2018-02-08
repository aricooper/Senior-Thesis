function [ G ] = geodesics(cs, u, handles)
%GEODESICS outputs 1X4 matrix of geodesic equations for each coordinate 
%takes christoffel symbols and four velocity to compute geodesics
G = sym(zeros(1,4));
%% Calculate acceleration

for k=1:4
    for i=1:4,
        for j=1:4;
            temp(j)  = -(cs(i,j,k)*u(i)*u(j));
        end
        temp1(i) = sym(sum(temp));
    end
    G(k) = simplify(sym(sum(temp1))); 
end

    
end

