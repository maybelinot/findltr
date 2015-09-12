import sys
import shelve


def group(young_lcp_parts, max_ltr_len, min_ltr_len, min_pattern_len, min_distance):
    '''FIXME: DOCS...'''
    groups_of_ltrs = [[[young_lcp_parts[0][0], young_lcp_parts[0][0] + min_pattern_len],
                       [young_lcp_parts[0][1], young_lcp_parts[0][1] + min_pattern_len]]]

    duplicates = False
    for lcp_part in young_lcp_parts[1:]:
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
    if duplicates or \
        (groups_of_ltrs[-1][0][1] - groups_of_ltrs[-1][0][0]) < min_ltr_len or \
            min_ltr_len > (groups_of_ltrs[-1][1][1] - groups_of_ltrs[-1][1][0]):
        del groups_of_ltrs[-1]
    return groups_of_ltrs
