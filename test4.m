clc;
clearvars;
% A small testcase for a recursive equation.
% The recursive equation is a_n = a_(n-1) * H_n / H_(n-1).
% H_n=psi(n+1)-psi(1) is the n-th term of the Harmonic sequence.
% The solution is obivously H_n itself.
% The guesses are x!, H_x, 1, x, ..., x^5 and pairwise products, such as 
% x^2*H_x.
syms x;
functions = {factorial(sym(x)), 2^sym(x), 3^sym(x), (psi(sym(x)+1)- psi(sym(1)))};
precision = 30; 
start = 0;
rec_degree = 3;
poly_degree = 7;
solution = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);

function ret = substitute(n, func_val)
    ret = (4*sym(n)^2-8*sym(n)+3)*func_val(:,1) - (5*sym(n)^2+1) * func_val(:,2) + (sym(n)^2-sym(n)+1) * func_val(:,3);
end
