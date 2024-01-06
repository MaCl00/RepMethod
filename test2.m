clc;
clearvars;

syms x;
functions = {gamma(sym(x)+1), psi(sym(x)+1)- psi(sym(1))};
precision = 300;
start = 0;
rec_degree = 3;
poly_degree = 3;
solution = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true); 


function ret = substitute(n, func_val)
    ret = func_val(:,1) - n * (n-1) * func_val(:,3);
end