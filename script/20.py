import gzip
import sys
import time
import argparse
from zbio import io, stat
from scipy.stats import chi2

parser = argparse.ArgumentParser()
#parser.add_argument('-i', dest='input', required=True, type=str, help='input vcf files')
#parser.add_argument('-o', dest='output', required=True, type=str, help='output vcf file')
parser.add_argument('-c', dest='chr', required=True, type=str, help='chr')
parser.add_argument('-s', dest='samples', type=str, default='alls_pass5.txt', help='Selected samples')
#parser.add_argument('-p', dest='pops', type=str, default='../allpops.txt', help='Sample populations')
parser.add_argument('--header', type=str, help='vcf headers to add')
parser.add_argument('--log', dest='log', type=str, help='output log dir')

args = parser.parse_args()

samples = []
for lst in io.splitIter(args.samples):
  samples.append(lst[0])

spop = {}
#for lst in io.splitIter(args.pops):
#  spop[lst[0]] = lst[2]
infos = {}
keptinfo = {'SVTYPE':1, 'SVLEN':1, 'END':1, 'MEI':1, '0START':1, '0END':1, 'MESTRAND':1, 'SOURCE':1, 'PERIOD':1, 'RU':1, }
for fn in ['sv2/1019_nyuwa_allsv_over05_integrate.vcf.gz']: # ['sv/1019_nyuwa_allsv_over05_integrate.vcf.gz']: # ('0830chuan/mea.vcf.gz','0830chuan/mei.vcf.gz', '0830chuan/SVchr{}.vcf.gz'.format(args.chr)):
  for lst in io.splitIter(fn):
    if lst[0].startswith('#'): continue
    if lst[6] not in ('PASS', '.'): continue
    chr = lst[0].replace('chr', '')
    if chr != args.chr: continue
    key = int(lst[1]), lst[2], lst[3], lst[4]
    out = []
    for s in lst[7].split(';'):
      lst2 = s.split('=')
      if lst2[0] in keptinfo: out.append(s)
    infos[key] = ';'.join(out)

trinfos = {}
for lst in io.splitIter('tr2/1019_nyuwa_tr_over05_or05_integrate.chr{}.vcf.gz'.format(args.chr)):
  if lst[0].startswith('#'): continue
  if lst[6] not in ('PASS', '.'): continue
  key = int(lst[1]), lst[2], lst[3]
  out = []
  for s in lst[7].split(';'):
    lst2 = s.split('=')
    if lst2[0] in keptinfo: out.append(s)
  trinfos[key] = ';'.join(out)

def flip_gt(s):
  lst = s.split('|')
  return '{}|{}'.format(lst[1], lst[0])

def chunhe_test(n, m, k, alt = 'two.tailed'): #pass
  p = 0
  for i in range(m+1):
    p += stat.hypergeo(n*2, n, m, i) * stat.hypergeo_test(n, i, m-i, k, alt=alt)
  return p

pv_rec = {}
def HWtest(AA, Aa, aa):
  n = float(AA + Aa + aa)
  if n == 0: return 1, 0
  pA = (2*AA + Aa)/(2*n)
  pa = 1 - pA # (Aa + 2*aa)/(2*n)
  EAA = n*(pA**2)
  EAa = n*(2*pA*pa)
  Eaa = n*(pa**2)
  if (EAA>=5) and (EAa>=5) and (Eaa>=5) and n>=40:
    chisqr = ((AA-EAA)**2)/EAA + ((Aa-EAa)**2)/EAa + ((aa-Eaa)**2)/Eaa
    p = chi2.sf(chisqr, 1)
  else:
    m = 2*aa + Aa
    n = int(n)
    if n not in pv_rec: pv_rec[n] = {}
    if (m, aa) not in pv_rec[n]:
      if aa > Eaa: p = 2 * chunhe_test(n, m, aa, alt='g')
      else: p = 2 * chunhe_test(n, m, aa, alt='l')
      if p > 1: p = 1
      pv_rec[n][(m, aa)] = p
    else: p = pv_rec[n][(m, aa)]
  return p, EAa

def gt2cnts(gts):
  data = [0, 0, 0]
  for gt in gts:
    if '.' in gt: continue
    elif gt in ('1|1', '1/1'): data[2] += 1
    elif gt in ('0|0', '0/0'): data[0] += 1
    elif gt in ('0|1', '1|0', '0/1', '1/0'): data[1] += 1
  return data

def vcf_iter(fn):
  mac, miss, gqth, hwth, hwth2 = 2, 0.9, 5, 1e-6, 1e-10
  n = c = 0
  for lst in io.splitIter(fn):
    if lst[0].startswith('##'): continue
    if lst[0].startswith('#'):
      ds = {}
      for i in range(9, len(lst)): ds[lst[i]] = i
      si = [ds[s] if s in ds else None for s in samples]
      if None in si:
        for s in samples:
          if s not in ds: sys.stderr.write('{} sample not in vcf {}\n'.format(fn, s))
    else:
      ost = lst[0:8] + ['GT']
      ost[0] = ost[0].replace('chr', '')
      if ost[0] != args.chr: continue
      #if lst[6] not in ('PASS', '.'):
      #  print('{} skip site not pass {}'.format(fn, lst[0:7]))
      #  continue
      ost[1] = int(ost[1])
      key = ost[1], ost[2], ost[3], ost[4]
      if key in infos: ost[7] = '{};{}'.format(infos[key], ost[7])
      else:
        key = ost[1], ost[2], ost[3]
        if key in trinfos: ost[7] = '{};{}'.format(trinfos[key], ost[7])
      if ost[3] == '.': ost[3] = 'N'
      n += 1
      gts = [lst[i].split(':')[0] if i is not None else './.' for i in si]
      #fmt = lst[8].split(':')
      #if False: # 'GQ' in fmt:
      #  qi = fmt.index('GQ')
      #  for i0, i in enumerate(si):
      #    gq = float(lst[i].split(':')[qi])
      #    if gq < gqth: gts[i0] = './.'
      l = len(si)
      #misth = l * miss
      nmis = ac = 0
      #data = gt2cnts(gts)
      #hw, ea = HWtest(data[0], data[1], data[2])
      #if ea < data[1] and hw < hwth:
      #  print('{} skip site HW ExHet {} {}|{} {}/{} {}'.format(fn, hw, data[1], round(ea, 3), c, n, lst[0:7]))
      #  continue
      c += 1
      yield ost + gts
  sys.stderr.write('{} end of vcf {}/{}\n'.format(fn, c, n))

#vis = [vcf_iter(fn) for fn in vcffiles]
#vs = [next(vi) for vi in vis]
infn = 'phasesv2/chr{}.sh2.vcf.gz'.format(args.chr)
vi = vcf_iter(infn)
v = next(vi)

#l = len(vcffiles)
#endi = [0] * l
l = len(samples)

skip = {}
'''
for lst in io.splitIter('merge/chr{}.same.txt'.format(args.chr)):
  if float(lst[0]) < 0.5: continue
  if lst[5] == '.': si = 0 # key = lst[3], lst[4], lst[5]
  elif lst[9] == '.': si = 1 #key = lst[7], lst[8], lst[9]
  elif lst[5].startswith('phase2'): si = 1
  elif lst[9].startswith('phase2'): si = 0
  else: continue
  key = lst[si*4+3], int(lst[si*4+4]), lst[si*4+5]
  skip[key] = 1
'''
header = ''
if args.header is not None:
  for line in open(args.header):
    header += line

#outfile = open(args.output, 'w')
outfile = sys.stdout
for lst in io.splitIter(infn):
  if lst[0].startswith('##'): outfile.write(io.tabjoin(lst, '\n'))
  elif lst[0].startswith('#'):
    outfile.write(header)
    outfile.write(io.tabjoin(lst[0:9], samples, '\n'))
    break


data, di = [], 0
#while sum(endi) < l:
for lst in io.splitIter('tr2phasing/split/chr{0}/tr.chr{0}.regions.txt'.format(args.chr)):
  sys.stderr.write(lst[0] + '\n')
  lst2 = lst[0].split('-')
  start, stop = int(lst2[0]), int(lst2[1])
  i = 0
  for i in range(len(data)):
    if data[i][1] >= start: break
  if i > 0:
    if di < i:
      for i1 in range(di, i): outfile.write(io.tabjoin(data[i1], '\n'))
      di = i
    data[0:i] = []
    di -= i

  while v is not None and v[1] <= stop:
    key = v[0], v[1], v[2]
    if key in skip:
      sys.stderr.write('Skip var {}\n'.format(v[0:5]))
    elif len(data) == 0 and v[1] < start: outfile.write(io.tabjoin(v, '\n'))
    else:
      data.append(v)
    try: v = next(vi)
    except StopIteration:
      v = None
      break
  tfn = 'tr2phasing/split/chr{0}/tr.chr{0}.{1}.phased.vcf.gz'.format(args.chr, lst[0])
  i = 0
  phase = [[0, 0] for j in range(l)]
  trdata = []
  for vt in vcf_iter(tfn):
    key = vt[1], vt[2], vt[3]
    if key in trinfos: # vt[2] == '.':
      trdata.append(vt)
      continue
    while i < len(data) and data[i][1] < vt[1]: i += 1
    while i < len(data) and data[i][1] == vt[1] and (data[i][3], data[i][4]) != (vt[3], vt[4]): i += 1
    #except:
    #  sys.stderr.write('Data {}, STR {}\n'.format(data[i][0:5], v[0:5]))
    #  exit(-1)
    if i >= len(data) or data[i][1] > vt[1]: continue
    for j in range(l):
      j1 = j + 9
      if vt[j1] in ('0|0', '1|1'): continue
      phase[j][0] += 1
      if data[i][j1] == vt[j1]: phase[j][1] += 1
      elif data[i][j1] == flip_gt(vt[j1]): phase[j][1] -= 1
      else: sys.stderr.write('Unknown GT {} {} {} {} {}\n'.format(j, vt[j1], data[i][j1], vt[0:5], data[i][0:5]))
  ps = ['{}|{}|{}'.format(j, phase[j][0], phase[j][1]) for j in range(l)]
  #sys.stderr.write(' '.join(ps) + '\n')
  for vt in trdata:
    while di < len(data) and data[di][1] <= vt[1]:
      outfile.write(io.tabjoin(data[di], '\n'))
      di += 1
    sys.stderr.write(io.tabjoin(vt[0:4], '\n'))
    ogt = [vt[j+9] if phase[j][1] >=0 else flip_gt(vt[j+9]) for j in range(l)]
    outfile.write(io.tabjoin(vt[0:9], ogt, '\n'))

for i1 in range(di, len(data)): outfile.write(io.tabjoin(data[i1], '\n'))
if v is not None:
  outfile.write(io.tabjoin(v, '\n'))
  for v in vi:
    outfile.write(io.tabjoin(v, '\n'))



