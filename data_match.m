function [match_index] = data_match(text, pattern, index, method)
% data_match - 数据匹配
%
% input:
%   - text: vector, 查找匹配数据, 长度为n
%   - pattern: vector, 待匹配数据, 长度为m, 需<n
%   - index: int, 匹配起始位置, 对于matlab, 从1开始
%   - method: str, 匹配方法, BF, BM, sunday, horspool等
% output:
%   - match_index: vector, 匹配数据的起始位置>0(matlab下标从1开始), 若为空, 则表示没有匹配成功
%
% usage: [match_index] = data_match(text, pattern, index, method)
%
% docs:
%   - 目前实现了上述几种方法
%

method = lower(method); % 转换成小写

% % switch方法
% switch method
% case 'bf'
%     [match_index] = BF(text, pattern, index);
% case 'bm'
%     [match_index] = BM(text, pattern, index);
% case 'sunday'
%     [match_index] = sunday(text, pattern, index);
% case 'horspool'
%     [match_index] = horspool(text, pattern, index);
% otherwise 
%     error(sprintf('%s is not correct parameter!', method));
% end

% 字典方法(matlab为结构体)
methods = struct('bf', @BF, 'bm', @BM, 'sunday', @sunday, 'horspool', @horspool);
if isfield(methods, method)
    func = getfield(methods, method);
    [match_index] = func(text, pattern, index);
else
    error(sprintf('%s is not correct parameter!', method));
end

end

%%
function [match_index] = BF(text, pattern, index)
% BF - 暴力匹配, 朴素算法
%
% input:
%   - text: vector, 查找匹配数据, 长度为n
%   - pattern: vector, 待匹配数据, 长度为m, 需<n
%   - index: int, 匹配起始位置, 对于matlab, 从1开始
% output:
%   - match_index: vector, 匹配数据的起始位置>0(matlab下标从1开始), 若为空, 则表示没有匹配成功
%
% docs：
%   - 暴力匹配, 就是穷举, 算法简单, 效率低
%   - pattern长度<text长度时，直接返回没有匹配
%   - 若有多个匹配, 则仅返回第一个匹配
%   - 对pattern, 在text中, 从头开始逐个匹配
%   - 最坏时间复杂度为O((n-m+1)*m)
%

% 默认值
match_index = [];

n = length(text);
m = length(pattern);

% 对匹配index进行限制
if index < 1
    index = 1;
end
if index > n
    index = n;
end

if n - index + 1 < m
    % 当数据长度没有pattern长时，直接返回
    return;
else
    % 逐个匹配
    i = index;
    while (i <= n-m+1)
        % 从左向右比较
        j = 1;
        while (j <= m && pattern(j) == text(i+j-1))
            j = j + 1;
        end
        if (j > m)
            % 找到匹配
            match_index = [match_index; i];
            offset = m;
        else
            % 没有找到匹配
            offset = 1;
        end
        i = i + offset;
    end
end

end

%%
function [match_index] = BM(text, pattern, index)
% BM - boyer moore算法(两个人的人名), 仅使用8bit数据
%
% input:
%   - text: vector, 查找匹配数据, 长度为n
%   - pattern: vector, 待匹配数据, 长度为m, 需<n
%   - index: int, 匹配起始位置, 对于matlab, 从1开始
% output:
%   - match_index: vector, 匹配数据的起始位置>0(matlab下标从1开始), 若为空, 则表示没有匹配成功
%
% docs：
%   - 从右向左匹配, 即从尾向头匹配
%   - 有2个匹配规则: 坏字符规则, 好后缀规则
%   - 坏字符: text中的数据与pattern不匹配时, 则text中的该数据为坏字符
%   - 好后缀: 在遇到坏字符之前, text和pattern已匹配成功的数据
%   - 坏字符规则: 当text中某个数据与pattern某个数据不匹配时, pattern需向右移动进行下一次匹配
%               移动的位数 = 坏字符在中的位置 - 坏字符在模式串中最右出现的位置, 若pattern中没有坏字符, 则最右出现位置为-1
%   - 好后缀规则: 当存在好后缀, 且失配时,
%               移动的位数 = 好后缀在pattern中对应的位置 - 好后缀在pattern中上一次出现的位置, 若好后缀没有出现在pattern中, 则上一次位置为-1
%   - 移动的位置 = max(坏字符, 好后缀)
%

% 默认值
match_index = [];

m = length(pattern);
n = length(text);

% 对匹配index进行限制
if index < 1
    index = 1;
end
if index > n
    index = n;
end

if n - index + 1 < m
    % 当数据长度没有pattern长时，直接返回
    return;
else
    bm_bc = bad_char(pattern);
    bm_gs = good_suffix(pattern);

    i = index;
    while (i <= n-m+1)
        % 从右向左比较, 找到坏字符位置
        j = m; % j为坏字符在pattern中的index
        while (j > 0 && pattern(j) == text(i+j-1))
            j = j - 1;
        end
        if j < 1
            % 找到匹配
            match_index = [match_index; i];
            offset = m;
            % offset = bm_gs(1);
        else
            % 没有找到匹配
            offset = max(bm_bc(text(i+j-1)) - m + j, bm_gs(j));
        end
        i = i + offset;
    end
end

end

function [bm_gs] = good_suffix(pattern)
% good_suffix - 计算好后缀(坏字符距离pattern尾部长度)
%
% input:
%   - pattern: vector, 待匹配数据, 长度为m
% output:
%   - bm_gs: vector, 坏字符距离pattern尾部的长度
%
% usage: [bm_gs] = good_suffix(pattern)
%
% docs:
%   - bm_gs(i)表示在pattern中index=i的数据(坏字符), 距离pattern尾部的长度
%   - case1: pattern中有子串与好后缀完全匹配
%   - case2: pattern中有子串与好后缀的字串匹配
%   - case3: pattern没有子串与好后缀及其字串匹配
%

m = length(pattern);

% 后缀长度数组
suff = suffix(pattern);

% case3
bm_gs = ones(m,1) * m;

% case2
j = 1;
for i = m:-1:1
    if (suff(i) == i)
        while j < m-i+1
            if (bm_gs(j) == m)
                bm_gs(j) = m - i;
            end
            j = j + 1;
        end
    end
end

% case1
for i = 1:m-1
    % suff(i)为好后缀的长度, 则对应的坏字符在pattern中的位置为m-suff(i)
    bm_gs(m-suff(i)) = m - i;
end

end

function [suff] = suffix(pattern)
% suffix - 公共后缀长度, 辅助计算
%
% input:
%   - pattern: vector, 待匹配数据, 长度为m
% output:
%   - suff: vector(256, 8bit), 后缀长度
%
% uaage: [suff] = suffix(pattern)
%
% docs:
%   - 处理时, 数据右对齐, 从右开始向左比较
%   - suff(i)表示以pattern(1:i)为后缀的数据, 与pattern(1:end)数据为后缀的公共后缀的长度, 即suff(i)=s, 使得pattern(i-s+1:i)=pattern(m-s+1:m)
%   - suff(i)>0时, 表示pattern有长度为suff(i)后缀在pattern(i)处可以找到相同的数据
%   - 如pattern = bcababab, 长度m=8
%   - 当 i=1 时, 以pattern(1:i)为后缀的数据为b, 以pattern(1:end)为后缀的数据为bcababab, 公共后缀为b, 长度为1
%   - 当 i=2 时, 以pattern(1:i)为后缀的数据为bc, 以pattern(1:end)为后缀的数据为bcababab, 公共后缀为b, 长度为0
%   ……
%   - 当 i=5 时, 以pattern(1:i)为后缀的数据为bcaba, 以pattern(1:end)为后缀的数据为bcababab, 没有公共后缀, 长度为0
%   - 当 i=6 时, 以pattern(1:i)为后缀的数据为bcabab, 以pattern(1:end)为后缀的数据为bcababab, 公共后缀为abab, 长度为4
%   - 当 i=7 时, 以pattern(1:i)为后缀的数据为bcababa, 以pattern(1:end)为后缀的数据为bcababab, 没有公共后缀, 长度为0
%   - 当 i=8 时, 以pattern(1:i)为后缀的数据为bcababab, 以pattern(1:end)为后缀的数据为bcababab, 公共后缀为bcababab, 长度为8
%

m = length(pattern);

suff = ones(m,1) * m;

% method1
% 从m-1开始
g = m;
for i = m-1:-1:1
    if (i > g && suff(i + m - f) < i - g)
        % g,f为上次匹配的index, 在[g,f]之间, 一定有pattern(i) = pattern(m-f+i)
        % 利用前面计算的suff
        suff(i) = suff(i + m - f);
    else
        if (i < g)
            g = i;
        end
        f = i; % f为pattern(1:i)右端的index, 为公共后缀尾部在pattern中的index
        while (g >= 1 && pattern(g) == pattern(g + m - f))
            % 从右向左比较, g为公共后缀头部在pattern中的index
            g = g - 1;
        end
        suff(i) = f - g;
    end
end

% % method2
% for i = m-1:-1:1
%     j = i;
%     while (j >= 1 && pattern(j) == pattern(j+m-i))
%         j = j - 1;
%     end
%     suff(i) = i - j;
% end

end

function [bm_bc] = bad_char(pattern)
% bad_char - 计算坏字符数组
%
% input:
%   - pattern: vector, 待匹配数据, 长度为m
% output:
%   - bm_bc: vector(长度为256, 8bit), 数据距离pattern尾部的长度
%
% usage: bm_bc = bad_char(pattern)
%
% docs:
%   - bm_bc(d)表示数据d距离pattern尾部的长度
%   - case1: 数据d没有出现在pattern中, bm_bc(d)=length(pattern)
%   - case2: 数据d出现在pattern中, bm_bc(d)=数据d在pattern中最后一次的位置距离pattern尾部的长度
%

m = length(pattern);
bm_bc = ones(256, 1) * m;

for index = 1 : m-1
    bm_bc(pattern(index)) = m - index;
end

end

%%
function [match_index] = sunday(text, pattern, index)
% sunday - sunday匹配算法
%
% input:
%   - text: vector, 查找匹配数据, 长度为n
%   - pattern: vector, 待匹配数据, 长度为m, 需<n
%   - index: int, 匹配起始位置, 对于matlab, 从1开始
% output:
%   - match_index: vector, 匹配数据的起始位置>0(matlab下标从1开始), 若为空, 则表示没有匹配成功
%
% docs:
%   - 从前向后(从左向右)匹配
%   - 当出现匹配失败时, 仅考虑text中未参与匹配的数据, 如text(i:i+m-1)与pattern进行匹配时, 考虑text(i+m), 有两种情况
%   - case1: text(i+m)不在pattern中, 移动位数=m+1
%   - case2: text(i+m)在pattern中, 移动位数=m-text(i+m)在pattern最右出现的index(以0开始)
%

% 默认值
match_index = [];

m = length(pattern);
n = length(text);

% 对匹配index进行限制
if index < 1
    index = 1;
end
if index > n
    index = n;
end

if n - index + 1 < m
    % 当数据长度没有pattern长时，直接返回
    return;
else
    offset_pre = sunday_pre(pattern);

    i = index;
    while (i <= n-m+1)
        % 从左向右比较
        j = 1;
        while (j <= m && pattern(j) == text(i+j-1))
            j = j + 1;
        end
        if (j > m)
            % 找到匹配
            match_index = [match_index; i];
            offset = m;
        else
            % 没有找到匹配
            if i+m > n
                % 已超出text范围
                return;
            end
            offset = offset_pre(text(i+m)+1);
        end
        i = i + offset;
    end
end

end

function [offset_pre] = sunday_pre(pattern)
% sunday_pre - 预处理, 不匹配时需要移动的位数, 仅8bit数据使用
%
% input:
%   - pattern: vector, 待匹配数据, 长度为m
% output:
%   - offset_pre: vector, 移动位数
%
% usage: [offset_pre] = sunday_pre(pattern)
%
% docs:
%   - 当出现匹配失败时, 仅考虑text中未参与匹配的数据, 如text(i:i+m-1)与pattern进行匹配时, 考虑text(i+m), 有两种情况
%   - case1: text(i+m)不在pattern中, 移动位数=m+1
%   - case2: text(i+m)在pattern中, 移动位数=m-text(i+m)在pattern最右出现的index(以0开始)
%
%

m = length(pattern);

% case1, 默认情况
offset_pre = ones(256,1) * (m+1);

% case2
for i = 1:m
    % matlab下表从1开始
    offset_pre(pattern(i) + 1) = m - i + 1;
end

end

%%
function [match_index] = horspool(text, pattern, index)
% horspool - horspool算法
%
% input:
%   - text: vector, 查找匹配数据, 长度为n
%   - pattern: vector, 待匹配数据, 长度为m, 需<n
%   - index: int, 匹配起始位置, 对于matlab, 从1开始
% output:
%   - match_index: vector, 匹配数据的起始位置>0(matlab下标从1开始), 若为空, 则表示没有匹配成功
%
% docs：
%   - 从右向左匹配, 即从尾向头匹配
%   - 当出现匹配失败时, 仅考虑text中与pattern对齐的最后一个数据, 如text(i:i+m-1)与pattern进行匹配时, 考虑text(i+m-1), 有两种情况
%   - case1: text(i+m-1)不在pattern前m-1个数据中, 移动位数=m
%   - case2: text(i+m-1)在pattern前m-1个数据中, 移动位数=m-1-text(i+m-1)在pattern最右出现的index(以0开始)
%

% 默认值
match_index = [];

m = length(pattern);
n = length(text);

% 对匹配index进行限制
if index < 1
    index = 1;
end
if index > n
    index = n;
end

if n - index + 1 < m
    % 当数据长度没有pattern长时，直接返回
    return;
else
    offset_pre = horspool_pre(pattern);

    i = index;
    while (i <= n-m+1)
        % 从右向左比较
        j = m;
        while (j > 0 && pattern(j) == text(i+j-1))
            j = j - 1;
        end
        if (j < 1)
            % 找到匹配
            match_index = [match_index; i];
            offset = m;
        else
            % 没有找到匹配
            if i+m > n
                % 已超出text范围
                return;
            end
            offset = offset_pre(text(i+m-1)+1);
        end
        i = i + offset;
    end
end

end

function [offset_pre] = horspool_pre(pattern)
% horspool_pre - 预处理, 不匹配时需要移动的位数, 仅8bit数据使用
%
% input:
%   - pattern: vector, 待匹配数据, 长度为m
% output:
%   - offset_pre: vector, 移动位数
%
% docs:
%   - 当出现匹配失败时, 仅考虑text中与pattern对齐的最后一个数据, 如text(i:i+m-1)与pattern进行匹配时, 考虑text(i+m-1), 有两种情况
%   - case1: text(i+m-1)不在pattern前m-1个数据中, 移动位数=m
%   - case2: text(i+m-1)在pattern前m-1个数据中, 移动位数=m-1-text(i+m-1)在pattern最右出现的index(以0开始)

m = length(pattern);

% case1, 默认情况
offset_pre = ones(256,1) * m;

% case2
for i = 1:m-1
    offset_pre(pattern(i) + 1) = m - i;
end

end

%%
