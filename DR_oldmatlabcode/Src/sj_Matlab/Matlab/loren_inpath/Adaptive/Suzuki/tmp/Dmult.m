function M = dmult(A,B)
%DMULT    DMULT(A,B) is the product of diag(A) and B

[mb,nb] = size(B);
M=(A*ones(1,nb)).*B;
