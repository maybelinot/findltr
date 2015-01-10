def original(seq, min_pattern_len, min_distance, max_distance):
    print('Original algorithm')
    output = []
    idx = 0
    while idx < len(seq) - (min_pattern_len * 2 + min_distance):
        pattern = seq[idx:idx + min_pattern_len]
        if not 'N' in pattern:
            text = seq[idx + min_pattern_len + min_distance:idx + min_pattern_len + min_distance + max_distance]
            if pattern in text:
                output.append([idx, text.index(pattern) + idx + min_pattern_len + min_distance])
                idx += min_pattern_len
        else:
            idx += min_pattern_len
        idx += 1

    return output

def kmp(seq, min_pattern_len, min_distance, max_distance):
    print('KnuthMorrisPratt algorithm')
    def KnuthMorrisPratt(text, pattern):
        pattern = list(pattern)
        shifts = [1] * (len(pattern) + 1)
        shift = 1
        for pos in range(len(pattern)):
            while shift <= pos and pattern[pos] != pattern[pos-shift]:
                shift += shifts[pos-shift]
            shifts[pos+1] = shift

        startPos = 0
        matchLen = 0
        for c in text:
            while matchLen == len(pattern) or \
                    matchLen >= 0 and pattern[matchLen] != c:
                startPos += shifts[matchLen]
                matchLen -= shifts[matchLen]
            matchLen += 1
            if matchLen == len(pattern):
                yield startPos

    output = []
    idx = 0
    while idx < len(seq) - (min_pattern_len * 2 + min_distance):
        pattern = seq[idx:idx + min_pattern_len]
        if not 'N' in pattern:
            text = seq[idx + min_pattern_len + min_distance:idx + min_pattern_len + min_distance + max_distance]
            ans = [el for el in KnuthMorrisPratt(text, pattern)]
            if ans:
                output.append([idx, ans[0] + idx + min_pattern_len + min_distance])
                idx += min_pattern_len
        else:
            idx += min_pattern_len
        idx += 1

    return output

def salcp(seq, min_pattern_len, min_distance, max_distance):
    print('Suffix array with long common prefix algorithm')
    def get_lcp(text, suffix_array):
        lcp = [0] * len(suffix_array)
        rank = [0] * len(suffix_array)
        n = len(suffix_array)
        for idx in range(n):
            rank[suffix_array[idx]] = idx
        h = 0
        for idx in range(n):
            k = rank[idx]
            if k == 0:
                lcp[k] = -1
            else:
                j = suffix_array[k - 1]
                while idx + h < n and j + h < n and text[idx + h] == text[j + h]:
                    h += 1
                lcp[k] = h
            if h > 0:
                h -= 1
        return lcp

    def count_lcp(text, pattern):
        min_len = min(len(text), len(pattern))
        if min_len < 1:
            return 0
        for i in range(min_len):
            if text[i] != pattern[i]:
                i -= 1
                break
        return i + 1

    def find_pattern(text, pattern, lcp, sa):
        left = 0
        right = len(text)
        middle = (right + left) // 2
        while True:
            if abs(right - left) <= 1:
                if pattern == text[sa[middle][1]:sa[middle][1] + len(pattern)]:
                    return sa[middle][1], len(pattern)
                return None
            pattern_lcp = max([count_lcp(pattern, text[sa[idx][1]:]) for idx in range(left, middle)])
            # lcp comparison
            if pattern_lcp > lcp[middle]:
                right = middle
                middle = (left + right) // 2
            elif pattern_lcp < lcp[middle]:
                left = middle
                middle = (left + right) // 2
            # pattern comparison
            elif len(pattern) <= len(text[sa[middle][1]:]):
                if pattern > text[sa[middle][1]:sa[middle][1] + len(pattern)]:
                    left = middle
                    middle = (left + right) // 2
                elif pattern < text[sa[middle][1]:sa[middle][1] + len(pattern)]:
                    right = middle
                    middle = (left + right) // 2
                else:
                    return sa[middle][1]
            else:
                if pattern[:len(text[sa[middle][1]:])] > text[sa[middle][1]:]:
                    left = middle
                    middle = (left + right) // 2
                elif pattern[:len(text[sa[middle][1]:])] < text[sa[middle][1]:]:
                    right = middle
                    middle = (left + right) // 2
                else:
                    return None

    output = []
    idx = 0
    while idx < len(seq) - (min_pattern_len * 2 + min_distance):
        pattern = seq[idx:idx + min_pattern_len]
        if not 'N' in pattern:
            text = seq[idx + min_pattern_len + min_distance:idx + min_pattern_len + min_distance + max_distance]
            suffix_array = sorted([(text[i:], i) for i in range(0, len(text))])
            lcp = get_lcp(text, list(map(lambda x: x[1], suffix_array)))
            ans = find_pattern(text, pattern, lcp, suffix_array)
            if None != ans:
                output.append([idx, ans+idx + min_pattern_len + min_distance])
                idx += min_pattern_len
        else:
            idx += min_pattern_len
        idx += 1
    
    return output