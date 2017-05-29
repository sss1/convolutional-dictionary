import numpy as np
import pyBigWig
import matplotlib.pyplot as plt

# bw = pyBigWig.open('data/wgEncodeBroadHistoneA549CtcfDex100nmSig.bigWig')
bw = pyBigWig.open('data/E123-H3K27ac.pval.signal.bigwig')
chromosome = 'chr21'
X = np.nan_to_num(np.asarray(bw.values(chromosome, 0, bw.chroms()[chromosome])))
# X = np.log(X)
# X[X == -np.inf] = 0
# plt.hist(X[np.logical_and(0.2 < X, X < 1)], 100)
# plt.show()

c = np.ceil(np.median(X))
n = 100
print 'Using cutoff ' + str(c) + ' and feature length ' + str(n)

to_keep = np.zeros(X.shape, dtype = bool)
last_nonzero_idx = -np.inf
for i in xrange(X.shape[0]): # traverse signal from left to right
  if X[i] >= c:
    last_nonzero_idx = i
  if i - last_nonzero_idx < n:
    to_keep[i] = True
next_nonzero_idx = np.inf
for i in reversed(xrange(X.shape[0])): # traverse signal from right to left
  if X[i] >= c:
    next_nonzero_idx = i
  if next_nonzero_idx - i < n:
    to_keep[i] = True

'''
to_keep now contains 1 at indices within n of an index above the cutoff, and
contains 0 everywhere else
'''

# f, axarr = plt.subplots(2, sharex=True)
# axarr[0].plot(X)
X = X[to_keep]
# axarr[1].plot(X)

# plt.plot(X)
# plt.show()

# for i in range(100):
#   c = float(i)/100
#   print str(float(np.count_nonzero(X > c)) / np.shape(X)[0]) + ' ' + str(c)

# TODO: Retain only samples within distance n of a sample above cutoff c.
# Use c = 0.19, since >98% of values are between 0.18 and 0.19.
np.savetxt("asCSV.csv", X, delimiter=",")
