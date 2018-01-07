load results.dat 
n = size(results,2)
A = results(1:n,:);
B = results(n+1:2*n,:);
C = results(2*n+1:3*n,:);
val = sum(sum(A*B-C,1),2);
if val == 0 
   disp("correct")
else
   disp("error")
endif
