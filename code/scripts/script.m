A=randn(3,3)

for i=2:3
    B(i-1,:)=cumsum(A(i,:)-A(i-1,:));
end