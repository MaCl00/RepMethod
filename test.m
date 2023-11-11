clc;
clearvars;
% A small testcase for a recursive equation.
% The recursive equation is a_n = a_(n-1) * H_n / H_(n-1).
% H_n=psi(n+1)-psi(1) is the n-th term of the Harmonic sequence.
% The solution is obivously H_n itself.
% The guesses are x!, H_x, 1, x, ..., x^5 and pairwise products, such as 
% x^2*H_x.
syms x;
functions = {gamma(sym(x)+1), psi(sym(x)+1)- psi(sym(1))};
precision = 100;
start = 0;
rec_degree = 2;
poly_degree = 5;
solution = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true); 


function ret = substitute(n, func_val)
    ret = func_val(1) - (psi(n+1)- psi(sym(1))) / (psi(n)- psi(sym(1))) * func_val(2);
end
