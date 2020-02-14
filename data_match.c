#include <stdio.h>
#include <string.h>

#define MAX_CHAR 256
#define SIZE 256
#define MAX(x, y) (x) > (y) ? (x) : (y)

void BF(char *pattern, int m, char *text, int n, int *match_index, int *mlen, int index);
void suffix(char *pattern, int m, int *suff);
void good_suffix(char *pattern, int m, int *bmgs);
void bad_char(char *pattern, int m, int *bmbc);
void BM(char *pattern, int m, char *text, int n, int *match_index, int *mlen, int index);
void sunday_pre(char *pattern, int m, int *offset_pre);
void sunday(char *pattern, int m, char *text, int n, int *match_index, int *mlen, int index);
void horspool_pre(char *pattern, int m, int *offset_pre);
void horspool(char *pattern, int m, char *text, int n, int *match_index, int *mlen, int index);
void print(int *array, int n, char *arrayName);

int main()
{
    char *text = "abcdefghijkcdeijk";
    char *pattern = "ijk";

    int mlen = 0;
    int match_index[10]; // 这里使用时需注意, 可能不止10个匹配

    // BF(pattern, strlen(pattern), text, strlen(text), match_index, &mlen, 0);
    // BM(pattern, strlen(pattern), text, strlen(text), match_index, &mlen, 0);
    // sunday(pattern, strlen(pattern), text, strlen(text), match_index, &mlen, 0);
    horspool(pattern, strlen(pattern), text, strlen(text), match_index, &mlen, 0);

    if (mlen)
    {
        print(match_index, mlen, "match_index:");
    }
    else
    {
        printf("no match!");
    }
}

void print(int *array, int n, char *arrayName)
{
    int i;
    printf("%s: ", arrayName);
    for (i = 0; i < n; i++)
    {
        printf("%d ", array[i]);
    }
    printf("\n");
}

/******************************************************************************
* Function   : BF, 暴力匹配算法
* Parameters :
*   - pattern       input   需匹配数据
*   - m             input   pattern数据长度
*   - text          input   源数据
*   - n             input   text数据长度
*   - match_index   output  匹配到的起始位置
*   - mlen          output  匹配到的个数, 如果为0表示没有匹配到
*   - index         input   开始匹配起始位置
* Returns    :
* Notes      :
* Changes    : 2020/02/12
******************************************************************************/
void BF(char *pattern, int m, char *text, int n, int *match_index, int *mlen, int index)
{
    // 对开始匹配起始位置进行限制
    if (index < 0)
    {
        index = 0;
    }
    if (index > n - 1)
    {
        index = n - 1;
    }

    if (n - index < m)
    {
        // 当数据长度没有pattern长时，直接返回
        return;
    }
    else
    {
        int i, j;
        int offset = 0;
        i = index;

        while (i <= n - m)
        {
            // 从左向右比较
            j = 0;
            while (j < m && text[i + j] == pattern[j])
            {
                j++;
            }
            if (j >= m)
            {
                // 找到匹配
                match_index[*mlen] = i;
                *mlen += 1;

                offset = m;
            }
            else
            {
                // 没有找到匹配
                if (i + m >= n)
                {
                    // 已超出text范围
                    return;
                }

                offset = 1;
            }
            i += offset; // 移动index
        }
    }
}

/******************************************************************************
* Function   : suffix
* Parameters :
*   - pattern   input   需匹配数据
*   - m         input   pattern数据长度
*   - suff      output  公共后缀的长度
* Returns    :
* Notes      :
*   - 处理时, 数据右对齐, 从右开始向左比较
*   - suff[i]表示以pattern[0:i]为后缀的数据, 与pattern[0:m-1]数据为后缀的公共后缀的长度, 即suff[i]=s, 使得pattern[i-s+1:i]=pattern[m-s+1:m]
*   - suff[i]>0时, 表示pattern有长度为suff[i]后缀在pattern(i)处可以找到相同的数据
*   - 如pattern = bcababab, 长度m=8
*   - 当 i=0 时, 以pattern[0:i]为后缀的数据为b, 以pattern[0:m-1]为后缀的数据为bcababab, 公共后缀为b, 长度为1
*   - 当 i=1 时, 以pattern[0:i]为后缀的数据为bc, 以pattern[0:m-1]为后缀的数据为bcababab, 公共后缀为b, 长度为0
*   ……
*   - 当 i=4 时, 以pattern[0:i]为后缀的数据为bcaba, 以pattern[0:m-1]为后缀的数据为bcababab, 没有公共后缀, 长度为0
*   - 当 i=5 时, 以pattern[0:i]为后缀的数据为bcabab, 以pattern[0:m-1]为后缀的数据为bcababab, 公共后缀为abab, 长度为4
*   - 当 i=6 时, 以pattern[0:i]为后缀的数据为bcababa, 以pattern[0:m-1]为后缀的数据为bcababab, 没有公共后缀, 长度为0
*   - 当 i=7 时, 以pattern[0:i]为后缀的数据为bcababab, 以pattern[0:m-1]为后缀的数据为bcababab, 公共后缀为bcababab, 长度为8
* Changes    : 2020/02/12
******************************************************************************/
void suffix(char *pattern, int m, int *suff)
{
    int f, g, i;

    suff[m - 1] = m;
    g = m - 1;

    for (i = m - 2; i >= 0; --i) // 从倒数第二个开始
    {
        if (i > g && suff[i + m - 1 - f] < i - g)
        {
            // g,f为上次匹配的index, 在[g,f]之间, 一定有pattern[i] = pattern[m-1-f+i]
            // 利用前面计算的suff
            suff[i] = suff[i + m - 1 - f];
        }
        else
        {
            if (i < g)
            {
                g = i;
            }
            f = i; // f为pattern[1:i]右端的index, 为公共后缀尾部在pattern中的index
            while (g >= 0 && pattern[g] == pattern[g + m - 1 - f])
            {
                // 从右向左比较, g为公共后缀头部在pattern中的index
                --g;
            }
            suff[i] = f - g;
        }
    }
}

/******************************************************************************
* Function   : good_suffix
* Parameters :
*   - pattern   input   需匹配数据
*   - m         input   pattern数据长度
*   - bmgs      output  坏字符距离pattern尾部的长度
* Returns    :
* Notes      :
*   - bmgs[i]表示在pattern中index=i的数据(坏字符), 距离pattern尾部的长度
*   - case1: pattern中有子串与好后缀完全匹配
*   - case2: pattern中有子串与好后缀的字串匹配
*   - case3: pattern没有子串与好后缀及其字串匹配
* Changes    : 2020/02/12
******************************************************************************/
void good_suffix(char *pattern, int m, int *bmgs)
{
    int i, j;

    int suff[SIZE];
    suffix(pattern, m, suff);

    // case3
    for (i = 0; i < m; i++)
    {
        bmgs[i] = m;
    }

    // case2
    j = 0;
    for (i = m - 1; i >= 0; i--)
    {
        if (suff[i] == i + 1)
        {
            for (; j < m - 1 - i; j++)
            {
                if (bmgs[j] == m)
                {
                    bmgs[j] = m - 1 - i;
                }
            }
        }
    }

    // case1
    for (i = 0; i <= m - 2; i++)
    {
        // suff[i]为好后缀的长度, 则对应的坏字符在pattern中的位置为m-1-suff[i]
        bmgs[m - 1 - suff[i]] = m - 1 - i;
    }
}

/******************************************************************************
* Function   : bad_char
* Parameters :
*   - pattern   input   需匹配数据
*   - m         input   pattern数据长度
*   - bmbc      output  坏字符距离pattern尾部的长度
* Returns    :
* Notes      :
*   - bmbc[d]表示数据d距离pattern尾部的长度
*   - case1: 数据d没有出现在pattern中, bmbc[d]=length(pattern)
*   - case2: 数据d出现在pattern中, bmbc[d]=数据d在pattern中最后一次的位置距离pattern尾部的长度
* Changes    : 2020/02/12
******************************************************************************/
void bad_char(char *pattern, int m, int *bmbc)
{
    int i;

    for (i = 0; i < MAX_CHAR; i++)
    {
        bmbc[i] = m;
    }

    for (i = 0; i < m - 1; i++)
    {
        bmbc[pattern[i]] = m - 1 - i;
    }
}

/******************************************************************************
* Function   : BM, boyer moore算法, 仅8bit数据使用
* Parameters :
*   - pattern       input   需匹配数据
*   - m             input   pattern数据长度
*   - text          input   源数据
*   - n             input   text数据长度
*   - match_index   output  匹配到的起始位置
*   - mlen          output  匹配到的个数, 如果为0表示没有匹配到
*   - index         input   开始匹配起始位置
* Returns    :
* Notes      :
*   - 从右向左匹配, 即从尾向头匹配
*   - 有2个匹配规则: 坏字符规则, 好后缀规则
*   - 坏字符: text中的数据与pattern不匹配时, 则text中的该数据为坏字符
*   - 好后缀: 在遇到坏字符之前, text和pattern已匹配成功的数据
*   - 坏字符规则: 当text中某个数据与pattern某个数据不匹配时, pattern需向右移动进行下一次匹配
*               移动的位数 = 坏字符在中的位置 - 坏字符在模式串中最右出现的位置, 若pattern中没有坏字符, 则最右出现位置为-1
*   - 好后缀规则: 当存在好后缀, 且失配时,
*               移动的位数 = 好后缀在pattern中对应的位置 - 好后缀在pattern中上一次出现的位置, 若好后缀没有出现在pattern中, 则上一次位置为-1
*   - 移动的位置 = max(坏字符, 好后缀)
* Changes    : 2020/02/12
******************************************************************************/
void BM(char *pattern, int m, char *text, int n, int *match_index, int *mlen, int index)
{
    // 对开始匹配起始位置进行限制
    if (index < 0)
    {
        index = 0;
    }
    if (index > n - 1)
    {
        index = n - 1;
    }

    if (n - index < m)
    {
        // 当数据长度没有pattern长时，直接返回
        return;
    }
    else
    {
        int i, j, bmbc[MAX_CHAR], bmgs[SIZE];
        int offset = 0;

        bad_char(pattern, m, bmbc);
        good_suffix(pattern, m, bmgs);

        i = index; // 从index开始匹配
        while (i <= n - m)
        {
            // 从右向左比较, 找到坏字符位置
            for (j = m - 1; j >= 0 && pattern[j] == text[i + j]; j--) // j为坏字符在pattern中的index
                ;
            if (j < 0)
            {
                // 找到匹配
                match_index[*mlen] = i;
                *mlen += 1;

                offset = bmgs[0];
            }
            else
            {
                // 没有找到匹配
                offset = MAX(bmbc[text[i + j]] - m + 1 + j, bmgs[j]);
            }

            i += offset; // 移动index
        }
    }
}

/******************************************************************************
* Function   : sunday_pre
* Parameters :
*   - pattern       input   需匹配数据
*   - m             input   pattern数据长度
*   - offset_pre    output  移动长度
* Returns    :
* Notes      :
*   - 当出现匹配失败时, 仅考虑text中未参与匹配的数据, 如text[i:i+m-1]与pattern进行匹配时, 考虑text[i+m], 有两种情况
*   - case1: text[i+m]不在pattern中, 移动位数=m+1
*   - case2: text[i+m]在pattern中, 移动位数=m-text[i+m]在pattern最右出现的index(以0开始)
* Changes    : 2020/02/12
******************************************************************************/
void sunday_pre(char *pattern, int m, int *offset_pre)
{
    int i;

    // case1, 默认
    for (i = 0; i < MAX_CHAR; i++)
    {
        offset_pre[i] = m + 1;
    }

    // case2
    for (i = 0; i < m; i++)
    {
        offset_pre[pattern[i]] = m - i;
    }
}

/******************************************************************************
* Function   : sunday
* Parameters :
*   - pattern       input   需匹配数据
*   - m             input   pattern数据长度
*   - text          input   源数据
*   - n             input   text数据长度
*   - match_index   output  匹配到的起始位置
*   - mlen          output  匹配到的个数, 如果为0表示没有匹配到
*   - index         input   开始匹配起始位置
* Returns    :
* Notes      :
*   - 从前向后(从左向右)匹配
*   - 当出现匹配失败时, 仅考虑text中未参与匹配的数据, 如text[i:i+m-1]与pattern进行匹配时, 考虑text[i+m], 有两种情况
*   - case1: text[i+m]不在pattern中, 移动位数=m+1
*   - case2: text[i+m]在pattern中, 移动位数=m-text[i+m]在pattern最右出现的index(以0开始)
* Changes    : 2020/02/12
******************************************************************************/
void sunday(char *pattern, int m, char *text, int n, int *match_index, int *mlen, int index)
{
    // 对开始匹配起始位置进行限制
    if (index < 0)
    {
        index = 0;
    }
    if (index > n - 1)
    {
        index = n - 1;
    }

    if (n - index < m)
    {
        // 当数据长度没有pattern长时，直接返回
        return;
    }
    else
    {
        int i, j, offset_pre[MAX_CHAR];
        int offset = 0;

        sunday_pre(pattern, m, offset_pre);

        i = index;
        while (i <= n - m)
        {
            // 从左向右比较
            j = 0;
            while (j < m && text[i + j] == pattern[j])
            {
                j++;
            }
            if (j >= m)
            {
                // 找到匹配
                match_index[*mlen] = i;
                *mlen += 1;

                offset = m;
            }
            else
            {
                // 没有找到匹配
                if (i + m >= n)
                {
                    // 已超出text范围
                    return;
                }

                offset = offset_pre[text[i + m]];
            }

            i += offset; // 移动index
        }
    }
}

/******************************************************************************
* Function   : horspool_pre
* Parameters :
*   - pattern       input   需匹配数据
*   - m             input   pattern数据长度
*   - offset_pre    output  移动长度
* Returns    :
* Notes      :
*   - 当出现匹配失败时, 仅考虑text中与pattern对齐的最后一个数据, 如text[i:i+m-1]与pattern进行匹配时, 考虑text[i+m-1], 有两种情况
*   - case1: text[i+m-1]不在pattern前m-1个数据中, 移动位数=m
*   - case2: text[i+m-1]在pattern前m-1个数据中, 移动位数=m-1-text[i+m-1]在pattern最右出现的index(以0开始)
* Changes    : 2020/02/12
******************************************************************************/
void horspool_pre(char *pattern, int m, int *offset_pre)
{
    int i;

    // case1, 默认
    for (i = 0; i < MAX_CHAR; i++)
    {
        offset_pre[i] = m;
    }

    // case2
    for (i = 0; i < m - 1; i++)
    {
        offset_pre[pattern[i]] = m - 1 - i;
    }
}

/******************************************************************************
* Function   : horspool
* Parameters :
*   - pattern       input   需匹配数据
*   - m             input   pattern数据长度
*   - text          input   源数据
*   - n             input   text数据长度
*   - match_index   output  匹配到的起始位置
*   - mlen          output  匹配到的个数, 如果为0表示没有匹配到
*   - index         input   开始匹配起始位置
* Returns    :
* Notes      :
*   - 从右向左匹配, 即从尾向头匹配
*   - 当出现匹配失败时, 仅考虑text中与pattern对齐的最后一个数据, 如text[i:i+m-1]与pattern进行匹配时, 考虑text[i+m-1], 有两种情况
*   - case1: text[i+m-1]不在pattern前m-1个数据中, 移动位数=m
*   - case2: text[i+m-1]在pattern前m-1个数据中, 移动位数=m-1-text[i+m-1]在pattern最右出现的index(以0开始)
* Changes    : 2020/02/12
******************************************************************************/
void horspool(char *pattern, int m, char *text, int n, int *match_index, int *mlen, int index)
{
    // 对开始匹配起始位置进行限制
    if (index < 0)
    {
        index = 0;
    }
    if (index > n - 1)
    {
        index = n - 1;
    }

    if (n - index < m)
    {
        // 当数据长度没有pattern长时，直接返回
        return;
    }
    else
    {
        int i, j, offset_pre[MAX_CHAR];
        int offset = 0;

        horspool_pre(pattern, m, offset_pre);

        i = index;
        while (i <= n - m)
        {
            // 从右向左比较
            j = m - 1;
            while (j >= 0 && text[i + j] == pattern[j])
            {
                j--;
            }
            if (j < 0)
            {
                // 找到匹配
                match_index[*mlen] = i;
                *mlen += 1;

                offset = m;
            }
            else
            {
                // 没有找到匹配
                if (i + m >= n)
                {
                    // 已超出text范围
                    return;
                }

                offset = offset_pre[text[i + m - 1]];
            }

            i += offset; // 移动index
        }
    }
}
