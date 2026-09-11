function values = sample_cast(values, options)
%SAMPLE_CAST Cast temporary ISDF sampling data without changing final solves.

if ~isfield(options, 'sample_precision') || ...
        strcmpi(options.sample_precision, 'double')
    return;
end
if strcmpi(options.sample_precision, 'single')
    values = single(values);
    return;
end
error('ISDF:SamplePrecision', ...
    'sample_precision must be ''double'' or ''single''.');
end
