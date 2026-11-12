from zbio import io
import sys
import numpy
import tabix

def gtc(s):
  s = s.split(':')[0]
  if '.' in s: return None
  return s.count('1')

def corr(d1, d2):
  if len(d1) != len(d2):
    print('Length not match! {} {}'.format(d1[0:5], d2[0:5]))
    exit(-1)
  m = {}
  gts = [[], []]
  for i in range(9, len(d1)):
    g1, g2 = gtc(d1[i]), gtc(d2[i])
    if g1 is None or g2 is None: continue
    gts[0].append(g1)
    gts[1].append(g2)
    key = g1, g2
    if key not in m: m[key] = 1
    else: m[key] += 1
  c = numpy.corrcoef(gts[0], gts[1])[0][1]
  return c, m

def getgts(lst):
  gts = []
  c, s = 0, 0
  for i in range(9, len(lst)):
    g = gtc(lst[i])
    gts.append(g)
    if g is None: continue
    c += 1
    s += g
  return gts, c, s

def corr2(gts1, gts2):
  l = len(gts1)
  gts = [[], []]
  for i in range(l):
    if i >= len(gts2): print(i, gts2)
    if gts1[i] is None or gts2[i] is None: continue
    gts[0].append(gts1[i])
    gts[1].append(gts2[i])
  c = numpy.corrcoef(gts[0], gts[1])[0][1]
  return c

#def gtcon(gts1, gts2):
#  l, c = len(gts1), 0.0
#  for i in range(l):
#    if gts1[i] is None or gts2[i] is None: continue
#    if gts1[i] == gts2[i]: c += 1
#    elif abs(gts1[i] - gts2[i]) == 1: c += 0.5
#  return c / l
def gtcon(gts1, gts2):
    matched = 0.0
    comp = 0
    for a, b in zip(gts1, gts2):
        if a is None or b is None: continue
        comp += 1
        if a == b:
            matched += 1.0
        elif abs(a - b) == 1:
            matched += 0.5
    return matched / comp if comp > 0 else float('nan')

#tmei = tabix.open('extract2/mei.vcf.gz')
#tmea = tabix.open('extract2/mea.vcf.gz')
t1k = tabix.open('/Parastor300s_G30S/haoxp/third/cmp/1019/extract/inv.mac1.967.vcf.gz') # extract2/63_sv_inv.vcf.gz')

outfile = open('/Parastor300s_G30S/haoxp/third/cmp/1019/svcmp.1kinv.txt', 'w')
outinv = open('/Parastor300s_G30S/haoxp/third/cmp/1019/svinfo.inv.txt', 'w')
#outdel = open('svinfo.del.txt', 'w')
out1k = open('/Parastor300s_G30S/haoxp/third/cmp/1019/svinfo.1kinv.txt', 'w')

ovlth, corth = 0.1, 0.5
for chr in range(1, 23):
  chr = str(chr)
  dins, ddel, dinv = [], [], []
  for lst in io.splitIter('/Parastor300s_G30S/haoxp/third/cmp/1019/extract/sv.mac1.chr{}.vcf.gz'.format(chr)):
    if lst[0].startswith('#'): continue
    if lst[4] in ('<DUP>', '<INS>'): data = dins
    elif lst[4] == '<DEL>': data = ddel
    elif lst[4] == '<INV>': data = dinv
    svlen = ''
    for s in lst[7].split(';'):
      if s.startswith('SVLEN='):
        svlen = s[6:]
        break
    if svlen == '.': continue
    if svlen == '': svlen = len(lst[4]) - len(lst[3])
    else: svlen = int(svlen)
    svlen = abs(svlen)
    pos = int(lst[1])
    gts, c, s = getgts(lst)
    data.append([pos, svlen, gts, c, s, lst[0:5]])

  data = dins
  for lst in []: # tmei.query('chr'+chr, 1, 300000000):
    if lst[6] != 'PASS': continue
    svlen = ''
    for s in lst[7].split(';'):
      if s.startswith('SVLEN='):
        svlen = s[6:]
        break
    if svlen == '.': continue
    svlen = int(svlen)
    pos = int(lst[1])
    gts, c, s = getgts(lst)
    data.append([pos, svlen, gts, c, s, lst[0:5]])

  data = ddel
  for lst in []: # tmea.query('chr'+chr, 1, 300000000):
    if lst[6] != 'PASS': continue
    svlen = ''
    for s in lst[7].split(';'):
      if s.startswith('SVLEN='):
        svlen = s[6:]
        break
    if svlen == '.': continue
    svlen = int(svlen)
    pos = int(lst[1])
    gts, c, s = getgts(lst)
    data.append([pos, svlen, gts, c, s, lst[0:5]])

  dins.sort(key=lambda di: di[0])
  ddel.sort(key=lambda di: di[0])
  dinv.sort(key=lambda di: di[0])
  
  for lst in t1k.query('chr'+chr, 1, 300000000): # io.splitIter('1019/final-vcf.unphased.vcf.gz'):
    if lst[0].startswith('#'): continue
    if lst[2].startswith('chr'): continue ### INS INV
  #if lst[2].startswith('phase2'):
  #  if lst[4].startswith('<INS'): ins = True
  #  else: ins = False
  #elif lst[4] in ('<DUP>', '<INS>'): ins = True
  #elif lst[4] == '<DEL>': ins = False
    l1, l2 = len(lst[3]), len(lst[4])
    data = dinv
    outd = outinv
    #if lst[2].find('INS') > 0:
    #  data = dins
    #  outd = outins
    #elif lst[2].find('DEL') > 0:
    #  data = ddel # ins = False
    #  outd = outdel
    #else:
    #  d = len(lst[4]) - len(lst[3])
    #  if d > 0:
    #    data = dins
    #    outd = outins
    #  else:
    #    data = ddel
    #    outd = outdel
  #elif len(lst[3]) >= 50: ins = False
  #elif len(lst[4]) >= 50: ins = True
  #else: continue
    #svlen = max(l1, l2) # ''
    end = ''
    for s in lst[7].split(';'):
      if s.startswith('END='):
        end = int(s[4:])
        break
    if end == '': continue
    #if svlen == '.': continue
    #if svlen == '': svlen = len(lst[4]) - len(lst[3])
  #else: svlen = int(svlen)
    #svlen = abs(int(svlen))
    pos = int(lst[1])
    svlen = end - pos + 1
    gts, c, s = getgts(lst)
  #if ins: data = dins
  #else: data = ddel
    i = 0
    for i in range(len(data)):
      if data[i][0] + data[i][1] < pos: continue
      break
    i0 = i
    for i in range(i0, len(data)):
      di = data[i]
      d = min(data[i][0] + data[i][1], end) - max(data[i][0], pos)
      if d <= 0: continue
      if d <= ovlth * data[i][1] or d <= ovlth * svlen: continue
      dr = min(d / data[i][1], d / svlen)
      if s == 0 or di[4] == 0: cor = 'nan'
      else: cor = corr2(data[i][2], gts)
      con = gtcon(data[i][2], gts)
      #if abs(c) >= corth:
      #  sins = 'Ins'
      # if not ins: sins = 'Del'
      outfile.write(io.tabjoin(lst[0:3], svlen, c, s, cor, con, data[i][5], di[1], di[3], di[4], dr, '\n'))
  #data[0:i0] = []
    for i in range(i0):
      di = data[i]
      outd.write(io.tabjoin(di[5], di[1], di[3], di[4], '\n'))
    data[0:i0] = []
    out1k.write(io.tabjoin(lst[0:3], svlen, c, s, '\n'))
  #data.append([pos, svlen, lst])
  for di in dinv: #, ddel, dinv
    outinv.write(io.tabjoin(di[5], di[1], di[3], di[4], '\n'))
  #for di in ddel:
  #  outdel.write(io.tabjoin(di[5], di[1], di[3], di[4], '\n'))


