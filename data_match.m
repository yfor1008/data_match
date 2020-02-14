function [match_index] = data_match(text, pattern, index, method)
% data_match - ����ƥ��
%
% input:
%   - text: vector, ����ƥ������, ����Ϊn
%   - pattern: vector, ��ƥ������, ����Ϊm, ��<n
%   - index: int, ƥ����ʼλ��, ����matlab, ��1��ʼ
%   - method: str, ƥ�䷽��, BF, BM, sunday, horspool��
% output:
%   - match_index: vector, ƥ�����ݵ���ʼλ��>0(matlab�±��1��ʼ), ��Ϊ��, ���ʾû��ƥ��ɹ�
%
% usage: [match_index] = data_match(text, pattern, index, method)
%
% docs:
%   - Ŀǰʵ�����������ַ���
%

method = lower(method); % ת����Сд

% % switch����
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

% �ֵ䷽��(matlabΪ�ṹ��)
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
% BF - ����ƥ��, �����㷨
%
% input:
%   - text: vector, ����ƥ������, ����Ϊn
%   - pattern: vector, ��ƥ������, ����Ϊm, ��<n
%   - index: int, ƥ����ʼλ��, ����matlab, ��1��ʼ
% output:
%   - match_index: vector, ƥ�����ݵ���ʼλ��>0(matlab�±��1��ʼ), ��Ϊ��, ���ʾû��ƥ��ɹ�
%
% docs��
%   - ����ƥ��, �������, �㷨��, Ч�ʵ�
%   - pattern����<text����ʱ��ֱ�ӷ���û��ƥ��
%   - ���ж��ƥ��, ������ص�һ��ƥ��
%   - ��pattern, ��text��, ��ͷ��ʼ���ƥ��
%   - �ʱ�临�Ӷ�ΪO((n-m+1)*m)
%

% Ĭ��ֵ
match_index = [];

n = length(text);
m = length(pattern);

% ��ƥ��index��������
if index < 1
    index = 1;
end
if index > n
    index = n;
end

if n - index + 1 < m
    % �����ݳ���û��pattern��ʱ��ֱ�ӷ���
    return;
else
    % ���ƥ��
    i = index;
    while (i <= n-m+1)
        % �������ұȽ�
        j = 1;
        while (j <= m && pattern(j) == text(i+j-1))
            j = j + 1;
        end
        if (j > m)
            % �ҵ�ƥ��
            match_index = [match_index; i];
            offset = m;
        else
            % û���ҵ�ƥ��
            offset = 1;
        end
        i = i + offset;
    end
end

end

%%
function [match_index] = BM(text, pattern, index)
% BM - boyer moore�㷨(�����˵�����), ��ʹ��8bit����
%
% input:
%   - text: vector, ����ƥ������, ����Ϊn
%   - pattern: vector, ��ƥ������, ����Ϊm, ��<n
%   - index: int, ƥ����ʼλ��, ����matlab, ��1��ʼ
% output:
%   - match_index: vector, ƥ�����ݵ���ʼλ��>0(matlab�±��1��ʼ), ��Ϊ��, ���ʾû��ƥ��ɹ�
%
% docs��
%   - ��������ƥ��, ����β��ͷƥ��
%   - ��2��ƥ�����: ���ַ�����, �ú�׺����
%   - ���ַ�: text�е�������pattern��ƥ��ʱ, ��text�еĸ�����Ϊ���ַ�
%   - �ú�׺: ���������ַ�֮ǰ, text��pattern��ƥ��ɹ�������
%   - ���ַ�����: ��text��ĳ��������patternĳ�����ݲ�ƥ��ʱ, pattern�������ƶ�������һ��ƥ��
%               �ƶ���λ�� = ���ַ����е�λ�� - ���ַ���ģʽ�������ҳ��ֵ�λ��, ��pattern��û�л��ַ�, �����ҳ���λ��Ϊ-1
%   - �ú�׺����: �����ںú�׺, ��ʧ��ʱ,
%               �ƶ���λ�� = �ú�׺��pattern�ж�Ӧ��λ�� - �ú�׺��pattern����һ�γ��ֵ�λ��, ���ú�׺û�г�����pattern��, ����һ��λ��Ϊ-1
%   - �ƶ���λ�� = max(���ַ�, �ú�׺)
%

% Ĭ��ֵ
match_index = [];

m = length(pattern);
n = length(text);

% ��ƥ��index��������
if index < 1
    index = 1;
end
if index > n
    index = n;
end

if n - index + 1 < m
    % �����ݳ���û��pattern��ʱ��ֱ�ӷ���
    return;
else
    bm_bc = bad_char(pattern);
    bm_gs = good_suffix(pattern);

    i = index;
    while (i <= n-m+1)
        % ��������Ƚ�, �ҵ����ַ�λ��
        j = m; % jΪ���ַ���pattern�е�index
        while (j > 0 && pattern(j) == text(i+j-1))
            j = j - 1;
        end
        if j < 1
            % �ҵ�ƥ��
            match_index = [match_index; i];
            offset = m;
            % offset = bm_gs(1);
        else
            % û���ҵ�ƥ��
            offset = max(bm_bc(text(i+j-1)) - m + j, bm_gs(j));
        end
        i = i + offset;
    end
end

end

function [bm_gs] = good_suffix(pattern)
% good_suffix - ����ú�׺(���ַ�����patternβ������)
%
% input:
%   - pattern: vector, ��ƥ������, ����Ϊm
% output:
%   - bm_gs: vector, ���ַ�����patternβ���ĳ���
%
% usage: [bm_gs] = good_suffix(pattern)
%
% docs:
%   - bm_gs(i)��ʾ��pattern��index=i������(���ַ�), ����patternβ���ĳ���
%   - case1: pattern�����Ӵ���ú�׺��ȫƥ��
%   - case2: pattern�����Ӵ���ú�׺���ִ�ƥ��
%   - case3: patternû���Ӵ���ú�׺�����ִ�ƥ��
%

m = length(pattern);

% ��׺��������
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
    % suff(i)Ϊ�ú�׺�ĳ���, ���Ӧ�Ļ��ַ���pattern�е�λ��Ϊm-suff(i)
    bm_gs(m-suff(i)) = m - i;
end

end

function [suff] = suffix(pattern)
% suffix - ������׺����, ��������
%
% input:
%   - pattern: vector, ��ƥ������, ����Ϊm
% output:
%   - suff: vector(256, 8bit), ��׺����
%
% uaage: [suff] = suffix(pattern)
%
% docs:
%   - ����ʱ, �����Ҷ���, ���ҿ�ʼ����Ƚ�
%   - suff(i)��ʾ��pattern(1:i)Ϊ��׺������, ��pattern(1:end)����Ϊ��׺�Ĺ�����׺�ĳ���, ��suff(i)=s, ʹ��pattern(i-s+1:i)=pattern(m-s+1:m)
%   - suff(i)>0ʱ, ��ʾpattern�г���Ϊsuff(i)��׺��pattern(i)�������ҵ���ͬ������
%   - ��pattern = bcababab, ����m=8
%   - �� i=1 ʱ, ��pattern(1:i)Ϊ��׺������Ϊb, ��pattern(1:end)Ϊ��׺������Ϊbcababab, ������׺Ϊb, ����Ϊ1
%   - �� i=2 ʱ, ��pattern(1:i)Ϊ��׺������Ϊbc, ��pattern(1:end)Ϊ��׺������Ϊbcababab, ������׺Ϊb, ����Ϊ0
%   ����
%   - �� i=5 ʱ, ��pattern(1:i)Ϊ��׺������Ϊbcaba, ��pattern(1:end)Ϊ��׺������Ϊbcababab, û�й�����׺, ����Ϊ0
%   - �� i=6 ʱ, ��pattern(1:i)Ϊ��׺������Ϊbcabab, ��pattern(1:end)Ϊ��׺������Ϊbcababab, ������׺Ϊabab, ����Ϊ4
%   - �� i=7 ʱ, ��pattern(1:i)Ϊ��׺������Ϊbcababa, ��pattern(1:end)Ϊ��׺������Ϊbcababab, û�й�����׺, ����Ϊ0
%   - �� i=8 ʱ, ��pattern(1:i)Ϊ��׺������Ϊbcababab, ��pattern(1:end)Ϊ��׺������Ϊbcababab, ������׺Ϊbcababab, ����Ϊ8
%

m = length(pattern);

suff = ones(m,1) * m;

% method1
% ��m-1��ʼ
g = m;
for i = m-1:-1:1
    if (i > g && suff(i + m - f) < i - g)
        % g,fΪ�ϴ�ƥ���index, ��[g,f]֮��, һ����pattern(i) = pattern(m-f+i)
        % ����ǰ������suff
        suff(i) = suff(i + m - f);
    else
        if (i < g)
            g = i;
        end
        f = i; % fΪpattern(1:i)�Ҷ˵�index, Ϊ������׺β����pattern�е�index
        while (g >= 1 && pattern(g) == pattern(g + m - f))
            % ��������Ƚ�, gΪ������׺ͷ����pattern�е�index
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
% bad_char - ���㻵�ַ�����
%
% input:
%   - pattern: vector, ��ƥ������, ����Ϊm
% output:
%   - bm_bc: vector(����Ϊ256, 8bit), ���ݾ���patternβ���ĳ���
%
% usage: bm_bc = bad_char(pattern)
%
% docs:
%   - bm_bc(d)��ʾ����d����patternβ���ĳ���
%   - case1: ����dû�г�����pattern��, bm_bc(d)=length(pattern)
%   - case2: ����d������pattern��, bm_bc(d)=����d��pattern�����һ�ε�λ�þ���patternβ���ĳ���
%

m = length(pattern);
bm_bc = ones(256, 1) * m;

for index = 1 : m-1
    bm_bc(pattern(index)) = m - index;
end

end

%%
function [match_index] = sunday(text, pattern, index)
% sunday - sundayƥ���㷨
%
% input:
%   - text: vector, ����ƥ������, ����Ϊn
%   - pattern: vector, ��ƥ������, ����Ϊm, ��<n
%   - index: int, ƥ����ʼλ��, ����matlab, ��1��ʼ
% output:
%   - match_index: vector, ƥ�����ݵ���ʼλ��>0(matlab�±��1��ʼ), ��Ϊ��, ���ʾû��ƥ��ɹ�
%
% docs:
%   - ��ǰ���(��������)ƥ��
%   - ������ƥ��ʧ��ʱ, ������text��δ����ƥ�������, ��text(i:i+m-1)��pattern����ƥ��ʱ, ����text(i+m), ���������
%   - case1: text(i+m)����pattern��, �ƶ�λ��=m+1
%   - case2: text(i+m)��pattern��, �ƶ�λ��=m-text(i+m)��pattern���ҳ��ֵ�index(��0��ʼ)
%

% Ĭ��ֵ
match_index = [];

m = length(pattern);
n = length(text);

% ��ƥ��index��������
if index < 1
    index = 1;
end
if index > n
    index = n;
end

if n - index + 1 < m
    % �����ݳ���û��pattern��ʱ��ֱ�ӷ���
    return;
else
    offset_pre = sunday_pre(pattern);

    i = index;
    while (i <= n-m+1)
        % �������ұȽ�
        j = 1;
        while (j <= m && pattern(j) == text(i+j-1))
            j = j + 1;
        end
        if (j > m)
            % �ҵ�ƥ��
            match_index = [match_index; i];
            offset = m;
        else
            % û���ҵ�ƥ��
            if i+m > n
                % �ѳ���text��Χ
                return;
            end
            offset = offset_pre(text(i+m)+1);
        end
        i = i + offset;
    end
end

end

function [offset_pre] = sunday_pre(pattern)
% sunday_pre - Ԥ����, ��ƥ��ʱ��Ҫ�ƶ���λ��, ��8bit����ʹ��
%
% input:
%   - pattern: vector, ��ƥ������, ����Ϊm
% output:
%   - offset_pre: vector, �ƶ�λ��
%
% usage: [offset_pre] = sunday_pre(pattern)
%
% docs:
%   - ������ƥ��ʧ��ʱ, ������text��δ����ƥ�������, ��text(i:i+m-1)��pattern����ƥ��ʱ, ����text(i+m), ���������
%   - case1: text(i+m)����pattern��, �ƶ�λ��=m+1
%   - case2: text(i+m)��pattern��, �ƶ�λ��=m-text(i+m)��pattern���ҳ��ֵ�index(��0��ʼ)
%
%

m = length(pattern);

% case1, Ĭ�����
offset_pre = ones(256,1) * (m+1);

% case2
for i = 1:m
    % matlab�±��1��ʼ
    offset_pre(pattern(i) + 1) = m - i + 1;
end

end

%%
function [match_index] = horspool(text, pattern, index)
% horspool - horspool�㷨
%
% input:
%   - text: vector, ����ƥ������, ����Ϊn
%   - pattern: vector, ��ƥ������, ����Ϊm, ��<n
%   - index: int, ƥ����ʼλ��, ����matlab, ��1��ʼ
% output:
%   - match_index: vector, ƥ�����ݵ���ʼλ��>0(matlab�±��1��ʼ), ��Ϊ��, ���ʾû��ƥ��ɹ�
%
% docs��
%   - ��������ƥ��, ����β��ͷƥ��
%   - ������ƥ��ʧ��ʱ, ������text����pattern��������һ������, ��text(i:i+m-1)��pattern����ƥ��ʱ, ����text(i+m-1), ���������
%   - case1: text(i+m-1)����patternǰm-1��������, �ƶ�λ��=m
%   - case2: text(i+m-1)��patternǰm-1��������, �ƶ�λ��=m-1-text(i+m-1)��pattern���ҳ��ֵ�index(��0��ʼ)
%

% Ĭ��ֵ
match_index = [];

m = length(pattern);
n = length(text);

% ��ƥ��index��������
if index < 1
    index = 1;
end
if index > n
    index = n;
end

if n - index + 1 < m
    % �����ݳ���û��pattern��ʱ��ֱ�ӷ���
    return;
else
    offset_pre = horspool_pre(pattern);

    i = index;
    while (i <= n-m+1)
        % ��������Ƚ�
        j = m;
        while (j > 0 && pattern(j) == text(i+j-1))
            j = j - 1;
        end
        if (j < 1)
            % �ҵ�ƥ��
            match_index = [match_index; i];
            offset = m;
        else
            % û���ҵ�ƥ��
            if i+m > n
                % �ѳ���text��Χ
                return;
            end
            offset = offset_pre(text(i+m-1)+1);
        end
        i = i + offset;
    end
end

end

function [offset_pre] = horspool_pre(pattern)
% horspool_pre - Ԥ����, ��ƥ��ʱ��Ҫ�ƶ���λ��, ��8bit����ʹ��
%
% input:
%   - pattern: vector, ��ƥ������, ����Ϊm
% output:
%   - offset_pre: vector, �ƶ�λ��
%
% docs:
%   - ������ƥ��ʧ��ʱ, ������text����pattern��������һ������, ��text(i:i+m-1)��pattern����ƥ��ʱ, ����text(i+m-1), ���������
%   - case1: text(i+m-1)����patternǰm-1��������, �ƶ�λ��=m
%   - case2: text(i+m-1)��patternǰm-1��������, �ƶ�λ��=m-1-text(i+m-1)��pattern���ҳ��ֵ�index(��0��ʼ)

m = length(pattern);

% case1, Ĭ�����
offset_pre = ones(256,1) * m;

% case2
for i = 1:m-1
    offset_pre(pattern(i) + 1) = m - i;
end

end

%%
