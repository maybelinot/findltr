__author__ = 'Эдуард'
dict = [{'text':'as','time':123},{'text':'ab','time':124},{'text':'bred','time':125},{'text':'berd','time':121},{'text':'aaaa','time':127}]
print(sorted(dict, key=lambda key:key['time']))