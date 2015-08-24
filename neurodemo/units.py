# Define a set of scaled unit symbols to make the code more clear
for unit in 'msVAΩFS':
    for pfx, val in [('p', -12), ('n', -9), ('u', -6), ('m', -3), ('c', -2), ('k', 3), ('M', 6), ('G', 9)]:
        locals()[pfx+unit] = 10**val
