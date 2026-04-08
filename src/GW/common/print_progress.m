function print_progress(current, total, varargin)
%PRINT_PROGRESS 显示当前时间和进度条，包含预计剩余时间
%   print_progress(current, total) 显示一个简单的进度条
%   print_progress(current, total, 'Message', msg) 显示带有自定义消息的进度条
%   print_progress(current, total, 'Message', msg, 'Task', taskName) 指定任务名称用于重置计时
%   print_progress(current, total, 'Total', totalAll, 'Message', msg) 在嵌套循环中使用总进度
%   print_progress(current, total, 'UpdateInterval', dt) 设置最小刷新间隔（秒，默认0.5s）
%   print_progress(current, total, 'PercentStep', p)   设置按百分比刷新（默认0，禁用）
%
%   输入参数：
%       current - 当前进度 (从1开始)
%       total   - 当前层级的总进度
%       'Message'       - 可选，进度条前的消息描述
%       'Task'          - 可选，任务标识符，用于区分不同任务的计时器
%       'Total'         - 可选，全局总进度（用于嵌套循环）
%       'Reset'         - 可选，true时重置计时器
%       'UpdateInterval'- 可选，最小刷新间隔（秒），默认0.5秒
%       'PercentStep'   - 可选，按百分比间隔刷新（如1表示每1%刷新一次），默认0（禁用）
%
%   示例：
%       % 每0.5秒刷新一次（默认）
%       print_progress(i, N, 'Message', 'Processing');
%       
%       % 每1秒刷新一次
%       print_progress(i, N, 'Message', 'Processing', 'UpdateInterval', 1.0);
%       
%       % 每1%进度刷新一次
%       print_progress(i, N, 'Message', 'Processing', 'PercentStep', 1);
%
%   示例输出：
%       [2026-04-01 17:12:00] Sigma           [████████████████░░░░░░░░░░░░░░░░░░░░]  40.00% | Elapsed: 00:02:15 | ETA: 00:03:23

persistent state;

% 解析可选参数
msg = '';
taskName = 'default';
totalAll = [];
resetTimer = false;
barWidth = 40; % 进度条宽度
minInterval = 0.5; % 默认最小刷新间隔（秒）
percentStep = 0;   % 默认不按百分比刷新

for i = 1:2:length(varargin)
    if length(varargin) < i+1, break; end
    param = varargin{i};
    if ischar(param) || isstring(param)
        switch lower(param)
            case 'message'
                msg = varargin{i+1};
            case 'task'
                taskName = varargin{i+1};
            case 'total'
                totalAll = varargin{i+1};
            case 'reset'
                resetTimer = true;
            case 'updateinterval'
                minInterval = varargin{i+1};
            case 'percentstep'
                percentStep = varargin{i+1};
        end
    end
end

% 判断是否是第一次调用新任务
isNewTask = isempty(state) || ~isfield(state, taskName) || ...
            isempty(state.(taskName)) || ~isfield(state.(taskName), 'startTime') || ...
            isempty(state.(taskName).startTime) || resetTimer;

% 初始化状态
if isNewTask
    if isempty(state)
        state = struct();
    end
    state.(taskName) = struct('startTime', tic(), 'lastUpdate', 0, 'lastPercent', 0, 'firstCallDone', false);
end

% 获取当前时间
currentTime = datetime('now');
% timeStr = datestr(currentTime, 'yyyy-mm-dd HH:MM:SS');
timeStr = datestr(currentTime, 'HH:MM:SS');

% 使用提供的总进度或当前层级的总进度
if ~isempty(totalAll)
    effectiveTotal = totalAll;
    effectiveCurrent = current;
else
    effectiveTotal = total;
    effectiveCurrent = current;
end

% 计算进度百分比
if effectiveTotal == 0
    percent = 0;
else
    percent = effectiveCurrent / effectiveTotal;
end

% 获取已用时间
elapsed = toc(state.(taskName).startTime);

% 判断是否需要更新
shouldUpdate = false;
if isNewTask || ~state.(taskName).firstCallDone
    % 第一次调用新任务，强制更新
    shouldUpdate = true;
    state.(taskName).firstCallDone = true;
else
    timeSinceLastUpdate = elapsed - state.(taskName).lastUpdate;
    percentSinceLastUpdate = (percent - state.(taskName).lastPercent) * 100;
    
    % 判断是否需要更新：时间间隔或百分比间隔
    if timeSinceLastUpdate >= minInterval
        shouldUpdate = true;
    end
    if percentStep > 0 && percentSinceLastUpdate >= percentStep
        shouldUpdate = true;
    end
    % 完成时总是更新
    if current >= total || effectiveCurrent >= effectiveTotal
        shouldUpdate = true;
    end
end

if ~shouldUpdate
    return;  % 跳过此次更新
end

state.(taskName).lastUpdate = elapsed;
state.(taskName).lastPercent = percent;

% 计算预计剩余时间 (ETA)
if effectiveCurrent > 1 && elapsed > 0.01 && percent > 0
    totalTime = elapsed / percent;
    remainingTime = totalTime - elapsed;
    etaStr = format_time(max(0, remainingTime));
else
    etaStr = '00:00:00';
end
elapsedStr = format_time(elapsed);

% 计算进度条
filledLength = round(barWidth * min(percent, 1));
emptyLength = barWidth - filledLength;

% 构建进度条 - 使用块状字符
if filledLength > 0
    bar = ['[', repmat('█', 1, filledLength), repmat('░', 1, emptyLength), ']'];
else
    bar = ['[', repmat('░', 1, emptyLength), ']'];
end

% 格式化输出
fprintf('\r'); % 回车符，覆盖当前行
if isempty(msg)
    fprintf('[%s] %s %6.2f%% | Elapsed: %s | ETA: %s', ...
        timeStr, bar, percent * 100, elapsedStr, etaStr);
else
    % 截断过长的消息以确保对齐
    maxMsgLen = 15;
    if length(msg) > maxMsgLen
        msg = msg(1:maxMsgLen);
    end
    fprintf('[%s] %-15s %s %6.2f%% | Elapsed: %s | ETA: %s', ...
        timeStr, msg, bar, percent * 100, elapsedStr, etaStr);
end

% 如果完成，换行
if current >= total || (effectiveCurrent >= effectiveTotal && effectiveCurrent > 0)
    fprintf('\n');
end

end

%% 辅助函数：格式化时间
function timeStr = format_time(seconds)
    if isempty(seconds) || seconds < 0
        seconds = 0;
    end
    if ~isfinite(seconds)
        seconds = 0;
    end
    hours = floor(seconds / 3600);
    minutes = floor(mod(seconds, 3600) / 60);
    secs = round(mod(seconds, 60));
    timeStr = sprintf('%02d:%02d:%02d', hours, minutes, secs);
end