function aqs = sigma_exchange_matrix_element(matrix_elements, nn, fallback)
%SIGMA_EXCHANGE_MATRIX_ELEMENT Return the bare-exchange element for band nn.

if isfield(matrix_elements, 'gme_exchange') && ...
        ~isempty(matrix_elements.gme_exchange)
    exchange = matrix_elements.gme_exchange;
    if isstruct(exchange)
        index = find(exchange.bands == nn, 1);
        if ~isempty(index)
            aqs = exchange.values(:, index);
            return;
        end
    else
        aqs = exchange(:, nn);
        return;
    end
end
aqs = fallback;
end
