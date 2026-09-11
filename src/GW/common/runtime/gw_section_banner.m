function gw_section_banner(section, phase)
%GW_SECTION_BANNER Print a consistent, readable GW calculation delimiter.

if nargin < 2
    phase = 'start';
end
section = upper(char(section));
switch lower(char(phase))
    case 'start'
        label = sprintf('BEGIN %s CALCULATION', section);
    case 'end'
        label = sprintf('END %s CALCULATION', section);
    otherwise
        error('gw_section_banner:InvalidPhase', ...
            'phase must be ''start'' or ''end''.');
end
line = repmat('=', 1, 78);
fprintf('\n%s\n%s\n%s\n', line, label, line);
end
