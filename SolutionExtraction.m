function solution = SolutionExtraction(matrix, x, precision, checkSolution)
    disp(size(matrix, 2))
    selectedColumns = [];
    remainingColumns = 1:size(matrix, 2);
    for i=1:size(matrix, 2)
        [~, max_index] = max(abs(x(remainingColumns)));
        most_influential_index = remainingColumns(max_index);
        remainingColumns(max_index) = [];
        selectedColumns = [selectedColumns, most_influential_index];
        reduced_matrix = matrix(:, selectedColumns);
        [~, ~, candidate_V] = svd(reduced_matrix);
        candidate_x = candidate_V(:, end);
        result = log10(norm(reduced_matrix * candidate_x, 'fro'));
        if result < -precision/2
            break;
        end
    end
    solution = zeros(size(matrix, 2), 1);
    solution(selectedColumns) = candidate_x;
    if (~isequal(checkSolution, 0) && size(candidate_x, 1) == size(matrix, 2)) || size(candidate_x, 1) == 1
        return;
    end
    solution(selectedColumns) = SolutionExtraction(reduced_matrix, candidate_x, precision, solution);
end