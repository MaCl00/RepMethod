clc;
clearvars;
syms x;
functions = {};
start = 0;
rec_degree = 3;
poly_degree = 5;
precision = 20;


for j=1:3
    homogeneous = 0;
    if j==1
        rec_eq = @substitute;
    elseif j==2
        rec_eq = @substitute2;
    else
        rec_eq = @substitute3;
        homogeneous = -(sym(x)-1)^4+4*(sym(x)-1)^2+5*sym(x)-3;
    end
    for k=1:3
        precision = 250*k;
        timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
        fid = fopen(['timings_', timestamp, '.txt'], 'w');
        for i=4:50
            [name, solution] = repertoire(functions, @substitute, 0, 10, start, rec_degree, 2, false, true);
            poly_degree = i; 
            disp([j,k,i]);
            tic
            [function_name, solution] = repertoire(functions, rec_eq, homogeneous, precision, start, rec_degree, poly_degree, false, true);
            elapsedTime = toc;
            fprintf(fid, '%d\t%f\n', precision, elapsedTime);
            disp(size(solution));
        end
        fclose(fid); % Close the file
    end
end

function ret = substitute(n, func_val)
    ret = sym(n)*(sym(n)-1)*func_val(:,1) - 2*(sym(n)+1)*(sym(n)-1)* func_val(:,2)+(sym(n)+1)*sym(n)*func_val(:, 3);
end
function ret = substitute2(n, func_val)
    ret = (sym(n)-2)*func_val(:,1) - (sym(n)-1)*func_val(:,2);
end
function ret = substitute3(n, func_val)
    ret = func_val(:,1) - (sym(n)-1)*func_val(:,2);
end
