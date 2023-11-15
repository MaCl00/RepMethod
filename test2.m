clc;
clearvars;

syms x;
functions = {gamma(sym(x)+1), psi(sym(x)+1)- psi(sym(1))};
precision = 200;
start = 0;
rec_degree = 2;
poly_degree = 3;
solution = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true); 


function ret = substitute(n, func_val)
    ret = func_val(:,1) - n * func_val(:,2);
end