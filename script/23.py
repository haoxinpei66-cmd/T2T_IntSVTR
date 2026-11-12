from zbio import io
import os
import sys
import argparse
import time
from os.path import isfile
#import tabix

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', required=True, type=str, help='input phased SNP vcf file')
parser.add_argument('-t', dest='str', required=True, type=str, help='input str vcf file')
parser.add_argument('-c', dest='chr', required=True, type=str, help='chr')
parser.add_argument('-o', dest='out', required=True, type=str, help='Out prefix')
#parser.add_argument('-s', dest='samples', type=str, help='Selected samples')
parser.add_argument('-f', dest='flank', type=int, default=50000, help='Flanking distance (default = 50 kb)')
#parser.add_argument('-m', dest='minaf', type=float, default=0.005, help='AF filter for SNPs (default = 0.005)')
#parser.add_argument('--hcdir', dest='hcdir', type=str, default='hapcut2', help='input hapcut2 dir')
#parser.add_argument('--use', type=str, help='Use only the listed variant sites')
parser.add_argument('--log', dest='log', type=str, help='output log dir')

args = parser.parse_args()
args.minafs = [0.05, 0.01, 0.005, 0.001, 0]

def unify_str(lst):
  print('unify_str', lst[1])
  ref = lst[3]
  infos = []
  anstr = ''
  for s in lst[7].split(';'):
    if s.startswith('AC='): acs = eval('[{}]'.format(s[3:]))
    elif s.startswith('AN='): an =  int(s[3:]); anstr = s
    else: infos.append(s)
  indata = {ref: [[0], an - sum(acs)]}
  i = 1
  for alt in lst[4].split(','):
    if alt not in indata: indata[alt] = [[i], acs[i-1]]
    else:
      indata[alt][0].append(i)
      indata[alt][1] += acs[i-1]
    i += 1
  d_by_len = {}
  for a in indata:
    if indata[a][1] == 0: continue
    l = len(a)
    if l not in d_by_len: d_by_len[l] = {}
    d_by_len[l][a] = indata[a][1]
  out, merge = {}, {}
  for l in d_by_len:
    d = d_by_len[l]
    if len(d) <= 1:
      for a in d: out[a] = d[a]
      continue
    dist = {}
    for a in d: dist[a] = {}
    dl = list(d.keys())
    for i in range(len(d)-1):
      for j in range(i+1, len(d)):
        m = 0
        for k in range(l):
          if dl[i][k] != dl[j][k]:
            m += 1
            if m > 1: break
        if m <= 1:
          dist[dl[i]][dl[j]] = m
          dist[dl[j]][dl[i]] = m
    change = False # True
    while change:
      change = False
      for a in dist:
        if a in merge: continue
        if a not in d: continue
        if len(dist[a]) != 1: continue # multiple links
        if a == ref: continue
        for a1 in dist[a]:
          if d[a] >= d[a1]: continue # not minor allele
          print('Merge {} {} into {} {}'.format(a, d[a], a1, d[a1]))
          d[a1] += d[a] # merge a into a1
          if a1 not in merge: merge[a1] = {}
          merge[a1][a] = d[a]
          del d[a]
          change = True
    for a in d: out[a] = d[a]
  ol = [[out[a], a] for a in out]
  ol.sort(reverse = True)
  #print(ol)
  oalts, oacs, aimap = [], [], {}
  i = 1
  for c, a in ol:
    #for i0 in indata[a][0]: aimap[str(i0)] = str(i)
    if a == ref: continue
    oalts.append(a)
    oacs.append(c)
    for i0 in indata[a][0]: aimap[str(i0)] = str(i)
    #if a in merge: aimap[str(i0)] = str(i)
    #for a1 in merge[a]: aimap[str(indata[a1][0])] = str(i)
    i += 1
  for i0 in indata[ref][0]: aimap[str(i0)] = '0'
  #if indata[ref][0] != 0: aimap[str(indata[ref][0])] = '0'
  #if ref in merge:
  #  for a1 in merge[ref]: aimap[str(indata[a1][0])] = '0'
  aimap['.'] = '.'
  print(aimap)

  rst = lst[:]
  rst[4] = ','.join(oalts)
  acstr = 'AC=' + ','.join(map(str, oacs))
  infos += [acstr, anstr]
  rst[7] = ';'.join(infos)
  oacs1 = [0] * len(oacs)
  oan = 0
  for i in range(9, len(lst)):
    if lst[i].startswith('.'): continue
    if lst[i].find('/') > 0: sep = '/'
    elif lst[i].find('|') > 0: sep = '|'
    else:
      print('GT error: {} {}'.format(i, lst[i]))
      exit(1)
    gts = lst[i].split(sep)
    ogts = aimap[gts[0]], aimap[gts[1]]
    rst[i] = '{}{}{}'.format(ogts[0], sep, ogts[1])
    if ogts[0] != '0': oacs1[int(ogts[0])-1] += 1
    if ogts[1] != '0': oacs1[int(ogts[1])-1] += 1
    oan += 2
  print('Input alleles: {}, output alleles: {}'.format(len(indata), len(out)))
  #print('Sum AC={}'.format(','.join(map(str, oacs))))
  print('Cnt AC={}, AN={}'.format(','.join(map(str, oacs1)), oan))
  return rst


def strIter(fn, isample, adjoin=5000):
  chr1 = 'chr' + args.chr
  l = len(isample)
  cluster, adj = [], adjoin
  imap = {}
  for lst in io.splitIter(fn):
    if lst[0].startswith('##'): continue
    if lst[0].startswith('#'):
      for i in range(9, len(lst)):
        s = lst[i]
        if s == 'Sample_WGC033697D': s = 'CN001381'
        if s not in isample: continue
        imap[isample[s]] = i
    if lst[0] != chr1: continue
    if lst[3] == '.' or lst[4] == '.': continue
    alts = lst[4].split(',')
    al = len(alts) + 1
    #isc = lst[8].split(':').index('SCORE')
    data = [0,0]
    acs = {i: 0 for i in range(al)}
    rst = lst[0:6] + ['PASS', lst[7], 'GT'] + ['./.'] * l
    c = 0
    for i in range(9, l+9):
      if i not in imap:
        #print('Sample not found in STR: {} {}'.format(i, sidx[i]))
        continue
      i2 = imap[i]
      if lst[i2].startswith('.'): continue
      lst2 = lst[i2].split(':')
      #gt, score = lst2[0].replace('|', '/'), float(lst2[isc])
      gt = lst2[0].replace('|', '/')
      #if score < 0.9: continue
      alleles = gt.split('/')
      a1, a2 = int(alleles[0]), int(alleles[1])
      acs[a1] += 1; acs[a2] += 1
      rst[i] = gt
      c += 1
    if c < l * 0.9: continue
    an, sumas, rare, hasz = 0, 0, True, False
    for i in range(al):
      an += acs[i]
      if acs[i] <= 2: sumas += acs[i]
      if acs[i] >= 2: rare = False
      if acs[i] == 0: hasz = True
    if rare: continue
    if sumas > l * 0.1 * 2: continue # too many singletons/doubletons

    rst[0] = rst[0].replace('chr', '')
    rst[1] = int(rst[1])
    if len(cluster) > 0 and rst[1] == cluster[-1][1]: continue ##
    #an = sum([acs[i] for i in range(al)])
    #rare = True
    #for i in range(1, al):
    #  if acs[i] >= 2:
    #    rare = False
    #    break
    #if rare: continue
    rst[7] += ';AC={};AN={}'.format(','.join([str(acs[i]) for i in range(1, al)]), an)
    #if hasz: 
    print(acs)
    rst = unify_str(rst)
    ralts = rst[4].split(',')
    ral = len(ralts) + 1
    if ral >= 50: # 120:
      print('Skipping multi-alleles >= 50 : {} sumas {}'.format(ral, sumas))
      #print('Merging multi-alleles > 120: {} sumas {}'.format(ral, sumas))
      print(io.tabjoin(rst[0:8]))
      continue
      #rst = unify_str(rst)
      #print(io.tabjoin(rst[0:8]))
      #ralts = rst[4].split(',')
      #ral = len(ralts) + 1
      #if ral > 126: continue
    if len(cluster) > 0 and cluster[-1][1] + adj < rst[1]:
      yield cluster
      cluster = []
    cluster.append(rst)
    adj = (args.flank - cluster[-1][1] + cluster[0][1]) * adjoin / args.flank
    if adj < 1000: adj = 1000
  yield cluster

isample, sidx = {}, {}
#outfile = sys.stdout
header = ''

for lst in io.splitIter(args.input):
  if lst[0].startswith('##'): header += io.tabjoin(lst, '\n') #outfile.write(io.tabjoin(lst, '\n'))
  elif lst[0].startswith('#'):
    header += io.tabjoin(lst, '\n')
    #outfile.write(io.tabjoin(lst, '\n'))
    for i, s in enumerate(lst):
      if i < 9: continue
      #if s == 'Sample_WGC033697D': s = 'CN001381'
      #if args.samples is not None and s not in args.samples: continue
      isample[s] = i
      sidx[i] = s
    break

regionfile = open('{}.chr{}.regions.txt'.format(args.out, args.chr), 'w')
c = 0
for cluster in strIter(args.str, isample):
  p1, p2 = cluster[0][1] - args.flank, cluster[-1][1] + args.flank
  if p1 < 0: p1 = 0
  outfn = '{}.chr{}.{}-{}.vcf'.format(args.out, args.chr, p1, p2)
  regionfile.write('{}-{}\n'.format(p1, p2))
  if isfile(outfn + '.gz'): continue ###
  #regionfile.write('{}-{}\n'.format(p1, p2))
  outfile = open(outfn, 'w')
  outfile.write(header)
  cmd = 'tabix {} {}:{}-{}'.format(args.input, args.chr, p1, p2)
  #print(cmd)
  res = os.popen(cmd)
  ti, tc = 0, len(cluster)
  #print('tc={}'.format(tc))
  n = 0
  minaflst = [[] for i in range(len(args.minafs))] # 1, ls2 = [], []
  for l in res.readlines():
    lst = l.split('\t', maxsplit=9)
    if len(lst) < 8: continue
    if len(lst[3]) > 50 or len(lst[4]) > 50: continue
    ac = an = 0
    lst2 = lst[7].split(';')
    for s in lst2:
      if s.startswith('AC='): ac = int(s.split('=')[1])
      elif s.startswith('AN='): an = int(s.split('=')[1])
    if an == 0: continue
    pos = int(lst[1])
    for i, minaf in enumerate(args.minafs):
      if ac >= an * minaf: minaflst[i].append([pos, l])# continue
    #if ac >= an * args.minaf2: ls2.append([pos, l])
  for i in range(len(args.minafs)):
    if len(minaflst[i]) >= 100: break
  #else: lines = ls2
  lines = minaflst[i]
  n = len(lines)
  for pos, l in lines:
    while ti < tc and pos > cluster[ti][1]:
      outfile.write(io.tabjoin(cluster[ti], '\n'))
      ti += 1
    outfile.write(l)
    #n += 1
  print('{}-{}, tc={}, n={}, afth={}'.format(p1, p2, tc, n, args.minafs[i]))
  while ti < tc:
    outfile.write(io.tabjoin(cluster[ti], '\n'))
    ti += 1
  outfile.close()
  os.popen('bgzip {0} && tabix {0}.gz'.format(outfn))
  c += 1
  #if c > 10:
  #  break

