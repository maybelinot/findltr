text = 'banana'
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

suffix_array = sorted([(text[i:], i) for i in range(0, len(text))])
lcp = get_lcp(text, list(map(lambda x: x[1], suffix_array)))
print(suffix_array)
print(lcp)