function repeated = epsilon_repeat_eden(eden, repeat_count)
%EPSILON_REPEAT_EDEN Match mapped gme column ordering.

if repeat_count == 1
    repeated = eden;
    return;
end
repeated = repmat(eden, repeat_count, 1, 1);
end
