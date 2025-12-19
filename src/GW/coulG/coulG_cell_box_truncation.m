function coulG = coulG_cell_box_truncation(g_q, gvec, sys)
%% set up fft_grid
Nfac = 3;
Nfft = zeros(1, 3);
Nfft(1) = gvec.n1;
Nfft(2) = gvec.n2;
Nfft(3) = gvec.n3;

for i = 1:3
    while ~check_FFT_size(Nfft(i), Nfac)
        Nfft(i) = Nfft(i) + 1;
    end
end

Nfft_double = 2 * Nfft;

for i = 1:3
    while ~check_FFT_size(Nfft_double(i), Nfac)
        Nfft_double(i) = Nfft_double(i) + 1;
    end
end

adot = inv(sys.bdot) * 4 * pi^2; % 实空间度规张量
for i = 1:3
    for j = 1:3
        adot(i, j) = adot(i, j) / (Nfft_double(i)*Nfft_double(j));
    end
end

scale = 2.0 * sqrt(det(adot));
trunc_shift = [0.5, 0.5, 0.5];
ncell = 3;
ws_dist = calculate_ws_distance(Nfft_double, adot, trunc_shift, ncell);
coul = scale ./ ws_dist;
coulG = fftn(coul);

%% From Nfft_double to Nfft
index_double = mod(g_q, Nfft_double) + 1;
linear_idx = sub2ind(Nfft_double, index_double(:, 1), index_double(:, 2), index_double(:, 3));
v_sampled = coulG(linear_idx);
Phase = 2 * pi * g_q ./ Nfft_double * trunc_shift';
coulG = real(v_sampled .* exp(-1i * Phase));

