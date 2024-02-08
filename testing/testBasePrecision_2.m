clc;
clearvars;
syms x;
functions = {};
starting = 2;
ending = 2000;
start = 0;
rec_degree = 3;
poly_degree = 5;
precision = 20;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
poly_degree = 10;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[~, ~] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
poly_degree = 15;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[~, ~] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
poly_degree = 5;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[~, ~] = repertoire(functions, @substitute2, -(sym(x)-1)^4+4*(sym(x)-1)^2+5*sym(x)-3, precision, start, rec_degree, poly_degree, false, true);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute2, -(sym(x)-1)^4+4*(sym(x)-1)^2+5*sym(x)-3, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
poly_degree = 10;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[~, ~] = repertoire(functions, @substitute2, -(sym(x)-1)^4+4*(sym(x)-1)^2+5*sym(x)-3, precision, start, rec_degree, poly_degree, false, true);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute2, -(sym(x)-1)^4+4*(sym(x)-1)^2+5*sym(x)-3, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
poly_degree = 15;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[~, ~] = repertoire(functions, @substitute2, -(sym(x)-1)^4+4*(sym(x)-1)^2+5*sym(x)-3, precision, start, rec_degree, poly_degree, false, true);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute2, -(sym(x)-1)^4+4*(sym(x)-1)^2+5*sym(x)-3, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
function ret = substitute(n, func_val)
    ret = (sym(n)-2)*func_val(:,1) - (sym(n)-1)*func_val(:,2);
end
function ret = substitute2(n, func_val)
    ret = func_val(:,1) - (sym(n)-1)*func_val(:,2);
end
