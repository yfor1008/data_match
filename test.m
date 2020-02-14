close all
clear
clc

text = 'abcdefghijkcdeijk';
pattern = 'cde';

text = abs(text);
pattern = abs(pattern);
% [match_index] = BF(text, pattern, 10);
% disp([match_index])
% [match_index] = BM(text, pattern, 1);
% disp([match_index])
% [match_index] = sunday(text, pattern, 1);
% disp([match_index])
% [match_index] = horspool(text, pattern, 1);
% disp([match_index])

[match_index] = data_match(text, pattern, 1, 'bm');
disp(match_index)
