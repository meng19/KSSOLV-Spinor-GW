function [c1, c2] = product_gram(left_components, right_components, ...
        ind_mu, idx_q, fftgrid, options)
%PRODUCT_GRAM Build product Gram factors without forming all products.
%   When IDX_Q and FFTGRID are supplied, C1 contains only the selected
%   Fourier samples of PP_mu'.  Its columns are formed and transformed in
%   bounded interpolation-point blocks instead of materializing ngrid*nmu.

if ~iscell(left_components)
    left_components = {left_components};
end
if ~iscell(right_components)
    right_components = {right_components};
end
if numel(left_components) ~= numel(right_components)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end
stream_to_g = nargin >= 5 && ~isempty(idx_q) && ~isempty(fftgrid);
if nargin < 6 || isempty(options)
    options = struct();
end

ngrid = size(left_components{1}, 1);
nleft = size(left_components{1}, 2);
nright = size(right_components{1}, 2);
for icomponent = 1:numel(left_components)
    if size(left_components{icomponent}, 1) ~= ngrid || ...
            size(right_components{icomponent}, 1) ~= ngrid
        error('ISDF:GridMismatch', ...
            'All left and right components must share the same grid size.');
    end
    if size(left_components{icomponent}, 2) ~= nleft || ...
            size(right_components{icomponent}, 2) ~= nright
        error('ISDF:BandMismatch', ...
            'All left and right components must share the same band counts.');
    end
end

ind_mu = ind_mu(:);
nmu = numel(ind_mu);
if stream_to_g
    idx_q = idx_q(:);
    c1 = complex(zeros(numel(idx_q), nmu, 'like', left_components{1}));
    block_size = max(1, floor(options.product_gram_block_elements / ...
        max(1, ngrid)));
else
    c1 = complex(zeros(ngrid, nmu, 'like', left_components{1}));
    block_size = nmu;
end
c2 = complex(zeros(nmu, nmu, 'like', left_components{1}));
left_mu = cell(size(left_components));
right_mu = cell(size(right_components));
for icomponent = 1:numel(left_components)
    left_mu{icomponent} = left_components{icomponent}(ind_mu, :);
    right_mu{icomponent} = right_components{icomponent}(ind_mu, :);
end
for icomponent = 1:numel(left_components)
    left_i_mu = left_mu{icomponent};
    right_i_mu = right_mu{icomponent};
    for jcomponent = 1:numel(left_components)
        left_j_mu = left_mu{jcomponent};
        right_j_mu = right_mu{jcomponent};
        c2 = c2 + (conj(left_i_mu) * left_j_mu.') .* ...
            (right_i_mu * conj(right_j_mu).');
    end
end

for first_mu = 1:block_size:nmu
    last_mu = min(nmu, first_mu + block_size - 1);
    mu_block = first_mu:last_mu;
    c1_block = complex(zeros(ngrid, numel(mu_block), ...
        'like', left_components{1}));
    for icomponent = 1:numel(left_components)
        left_i = left_components{icomponent};
        right_i = right_components{icomponent};
        for jcomponent = 1:numel(left_components)
            c1_block = c1_block + (conj(left_i) * ...
                left_mu{jcomponent}(mu_block, :).') .* ...
                (right_i * conj(right_mu{jcomponent}(mu_block, :)).');
        end
    end
    if stream_to_g
        c1(:, mu_block) = local_fft_samples(c1_block, idx_q, fftgrid);
    else
        c1(:, mu_block) = c1_block;
    end
end
end

function values_g = local_fft_samples(values_real, idx_q, fftgrid)
values_g = complex(zeros(numel(idx_q), size(values_real, 2), ...
    'like', values_real));
for icolumn = 1:size(values_real, 2)
    values_fft = ifftn(reshape(values_real(:, icolumn), fftgrid));
    values_g(:, icolumn) = values_fft(idx_q);
end
end
