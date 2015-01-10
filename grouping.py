import sys, shelve

print('Grouping of patterns')
max_ltr_len, min_ltr_len, min_pattern_len, min_distance = [int(item) for item in sys.argv[1:]]

db = shelve.open('LCP.db', writeback=True)
db['young_lcp'] = []
if db['young_lcp_parts']==[]:
        sys.exit(2)
groups_of_ltrs = [[[db['young_lcp_parts'][0][0], db['young_lcp_parts'][0][0] + min_pattern_len],
                   [db['young_lcp_parts'][0][1], db['young_lcp_parts'][0][1] + min_pattern_len]]]

duplicates = False
for lcp_part in db['young_lcp_parts'][1:]:
    if lcp_part[0] + min_pattern_len + min_distance > groups_of_ltrs[-1][1][0]:
        if lcp_part[0] > groups_of_ltrs[-1][1][1]:
            if duplicates\
                    or (groups_of_ltrs[-1][0][1] - groups_of_ltrs[-1][0][0]) < min_ltr_len \
                    or (groups_of_ltrs[-1][1][1] - groups_of_ltrs[-1][1][0]) < min_ltr_len:
                groups_of_ltrs[-1] = [[lcp_part[0], lcp_part[0] + min_pattern_len],
                                      [lcp_part[1], lcp_part[1] + min_pattern_len]]
                duplicates = False
            else:
                groups_of_ltrs.append([[lcp_part[0], lcp_part[0] + min_pattern_len],
                                       [lcp_part[1], lcp_part[1] + min_pattern_len]])
        else:
            duplicates = True
    elif (lcp_part[0] - groups_of_ltrs[-1][0][0] < max_ltr_len) or \
            (lcp_part[1] - groups_of_ltrs[-1][1][0] < max_ltr_len):
        groups_of_ltrs[-1][0][1] = lcp_part[0] + min_pattern_len
        groups_of_ltrs[-1][1][1] = lcp_part[1] + min_pattern_len
    else:
        duplicates = True
if duplicates or (groups_of_ltrs[-1][0][1] - groups_of_ltrs[-1][0][0]) < min_ltr_len \
        or min_ltr_len > (groups_of_ltrs[-1][1][1] - groups_of_ltrs[-1][1][0]):
    del groups_of_ltrs[-1]
db['young_lcp'] = groups_of_ltrs
db.close()