function [M, function_name] = generateGuesses(func_guess, range, num_poly, verbose)
    syms x;
    num_func = length(func_guess);
    if verbose
        fprintf('Evaluating function-guesses');
    end
    sequences = cellfun(@(j) cellfun(@(i) subs(func_guess{j}, x, sym(i)), range), num2cell(1:num_func), 'UniformOutput', false);
    if verbose
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        disp('Finished evaluating function-guesses');
        fprintf('Calculating product-functions');
    end
    M(1:num_func, :) = vertcat(sequences{:});
    current_row = num_func + 1;
    N = length(sequences);
    % Adding f*g if f and g are in the input
    counter = 1;
    prod_func = cell(1, (N*N-N)/2);
    for i = 1:N
        for j = i+1:N
            M(current_row, :) = M(i, :) .* M(j, :);
            current_row = current_row + 1;
            prod_func{counter} = func_guess{i} * func_guess{j};
            counter = counter + 1;
        end
    end
    func_guess = [func_guess, prod_func];
    if verbose
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        disp('Finished calculating product-functions');
        fprintf("Calculating polynomes");
    end
    % Calculating polynomes up to the degree "rec_degree"
    polynomes = cell(1, num_poly);
    poly_name = cell(1, num_poly);
    if num_poly >= 1
        polynomes{1} = cellfun(@(i) 1, range);
        poly_name{1} = 1;
    end
    for j=1:num_poly-1
        polynomes{j+1} = cellfun(@(i) polynomes{j}(i- range{1}(1) + 1) * (sym(i) - j + 1 - range{1}(1)) / j, range);
        poly_name{j+1} = prod(sym(x)-j+1- range{1}(1):sym(x)- range{1}(1)) / factorial(j);
    end
    polyMatrix = vertcat(polynomes{:});
    if verbose
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        disp('Finished calculating polynomes');
        fprintf("Calculating function-polynom products");
    end
    % Adding f*p for f is in the input and p is a polynom previous
    % calculated
    counter = 1;
    function_name = cell(1, (current_row-1)*(num_poly-1));
    for i = 1:(current_row-1)
        for j= 2:num_poly
            M(current_row, :) = M(i, :) .* polyMatrix(j, :);
            current_row = current_row + 1;
            function_name{counter} = func_guess{i} *poly_name{j};
            counter = counter + 1;
        end
    end
    M(current_row:current_row + num_poly - 1, :) = polyMatrix;
    function_name = [func_guess, function_name, poly_name];
    if verbose
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        disp("Finished calulating function-values");
    end
end