import itertools
from operator import mul

mutations = [64, 64, 32, 32, 128, 64]
library_size = 0

for combi in itertools.combinations(mutations, 3):
	combi_size = reduce(mul, combi, 1) # The product of the combi list
	library_size += combi_size
	print str(combi) + '\t' + str(combi_size)
	
print library_size