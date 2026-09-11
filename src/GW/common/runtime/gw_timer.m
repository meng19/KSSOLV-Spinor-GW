function varargout = gw_timer(action, name)
%GW_TIMER Lightweight CPU/WALL profiler for GW workflow sections.
%   gw_timer('reset'), gw_timer('start', NAME), gw_timer('stop', NAME),
%   gw_timer('report', TITLE).  Timers with the same NAME are accumulated.

persistent active stats
if isempty(active)
    active = struct();
    stats = struct();
end
switch lower(action)
    case 'reset'
        active = struct();
        stats = struct();
    case 'start'
        key = local_key(name);
        active.(key) = struct('wall', tic, 'cpu', cputime);
    case 'stop'
        key = local_key(name);
        if ~isfield(active, key)
            return;
        end
        wall = toc(active.(key).wall);
        cpu = cputime - active.(key).cpu;
        active = rmfield(active, key);
        if ~isfield(stats, key)
            stats.(key) = struct('name', name, 'calls', 0, ...
                'cpu_total', 0, 'cpu_min', inf, 'cpu_max', 0, ...
                'wall_total', 0, 'wall_min', inf, 'wall_max', 0);
        end
        entry = stats.(key);
        entry.calls = entry.calls + 1;
        entry.cpu_total = entry.cpu_total + cpu;
        entry.cpu_min = min(entry.cpu_min, cpu);
        entry.cpu_max = max(entry.cpu_max, cpu);
        entry.wall_total = entry.wall_total + wall;
        entry.wall_min = min(entry.wall_min, wall);
        entry.wall_max = max(entry.wall_max, wall);
        stats.(key) = entry;
    case 'report'
        if nargin < 2 || isempty(name)
            name = 'Timing information';
        end
        fprintf('\n%s\n\n', name);
        fprintf('%-28s %9s %9s %9s %9s %9s\n', ...
            'Routine', 'CPU min', 'CPU max', 'WALL min', 'WALL max', 'Calls');
        fprintf('%s\n', repmat('-', 1, 82));
        entries = struct2cell(stats);
        for i = 1:numel(entries)
            entry = entries{i};
            fprintf('%-28s %9.2f %9.2f %9.2f %9.2f %9d\n', ...
                entry.name, entry.cpu_min, entry.cpu_max, ...
                entry.wall_min, entry.wall_max, entry.calls);
        end
    otherwise
        error('GW:TimerAction', 'Unknown timer action "%s".', action);
end
end

function key = local_key(name)
key = matlab.lang.makeValidName(char(name));
end
