clc;
clearvars;
syms x;
functions = {(-1)^sym(x)};
precision = 5; 
start = 0;
rec_degree = 2;
poly_degree = 5;
[function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true, true);
function ret = substitute(n, func_val)
    ret = func_val(:,1) + sym((n+2)^4)/sym((n+1)^4)*func_val(:,2);
end
