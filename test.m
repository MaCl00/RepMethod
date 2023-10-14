clc;
clearvars;
syms x;
functions = {gamma(sym(x)+1), psi(sym(x)+1)- psi(sym(1))};
precision = 50;
start = 0;
rec_degree = 3;
poly_degree = 5;
solution = repertoire(functions, @substitute, precision, start, rec_degree, poly_degree, false); 

function ret = substitute(n, func_val)
    ret = (func_val(1) - (psi(n+1)- psi(sym(1))) / (psi(n)- psi(sym(1))) * func_val(2));
end
