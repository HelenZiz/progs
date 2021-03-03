import scipy.spatial.distance as sp
import re
n = []
d = dict()
k = 0
l = 1
ind = 0
file = open('dataset.txt')
for line in file:
    n.append(line.lower())
for i in range(len(n)):
    for word in re.split('[^a-z]', n[i]):
        if (word != '') & (word != ',') & (word != '.'):
            if word not in d:
                d[word] = k
                k += 1
m = [[0]*k for i in range(len(n))]
for j in range(len(n)):
    for word in re.split('[^a-z]', n[j]):
        if word in d:
            m[j][d[word]] += 1
s = [0]*(len(n)-1)
for j in range(len(n)-1):
    s[j] = sp.cosine(m[0], m[j+1])
for i in range(len(s)):
    if s[i]<s[ind]:
        l = ind
        ind = i
print(ind+1, l+1)
file.close()
