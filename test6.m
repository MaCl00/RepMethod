clc;
clearvars;
% A small testcase for a recursive equation.
% The recursive equation is a_n = a_(n-1) * H_n / H_(n-1).
% H_n=psi(n+1)-psi(1) is the n-th term of the Harmonic sequence.
% The solution is obivously H_n itself.
% The guesses are x!, H_x, 1, x, ..., x^5 and pairwise products, such as 
% x^2*H_x.
syms x;
functions = {factorial(sym(x)), (2)^sym(x), (-1)^sym(x), (psi(sym(x)+1)- psi(sym(1)))};
precision = 100; 
start = 0;
rec_degree = -1;
poly_degree = 1;
[function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true, true);
disp(solution)
function ret = substitute(~, func_val)
    ret = func_val(:,1) - sum(func_val(:,2:end), 2);
end
