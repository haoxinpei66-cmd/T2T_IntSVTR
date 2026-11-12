import gzip
import sys
import time
import argparse
from zbio import io, stat
from scipy.stats import chi2

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', required=True, type=str, help='input vcf files')
parser.add_argument('-o', dest='output', required=True, type=str, help='output vcf file')
parser.add_argument('-c', dest='chr', required=True, type=str, help='chr')
parser.add_argument('-s', dest='samples', type=str, default='alls_pass5.txt', help='Selected samples')
parser.add_argument('-p', dest='pops', type=str, default='../allpops.txt', help='Sample populations')
parser.add_argument('--header', type=str, help='vcf headers to add')
parser.add_argument('-f', dest='filt', type=str, help='Filter unpassed variant sites')
parser.add_argument('--use', type=str, help='Use only the listed variant sites')
parser.add_argument('--log', dest='log', type=str, help='output log dir')

args = parser.parse_args()
vcffiles = args.input.split(',')

samples = []
for lst in io.splitIter(args.samples):
  samples.append(lst[0])

spop = {}
for lst in io.splitIter(args.pops):
  spop[lst[0]] = lst[2]

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
          if s not in ds: print('{} sample not in vcf {}'.format(fn, s))
    else:
      ost = lst[0:8] + ['GT']
      ost[0] = ost[0].replace('chr', '')
      if ost[0] != args.chr: continue
      if lst[6] not in ('PASS', '.'):
        print('{} skip site not pass {}'.format(fn, lst[0:7]))
        continue
      ost[1] = int(ost[1])
      n += 1
      gts = [lst[i].split(':')[0] if i is not None else './.' for i in si]
      fmt = lst[8].split(':')
      if False: # 'GQ' in fmt:
        qi = fmt.index('GQ')
        for i0, i in enumerate(si):
          gq = float(lst[i].split(':')[qi])
          if gq < gqth: gts[i0] = './.'
      l = len(si)
      misth = l * miss
      nmis = ac = 0
      for gt in gts:
        if gt == './.': #'.' in gt:
          nmis += 1
          if nmis > misth: break
        elif '1' in gt: ac += 1
      if nmis > misth:
        print('{} skip site high miss rate >= {} {}/{} {}'.format(fn, nmis, c, n, lst[0:7]))
        continue
      if ac < mac:
        print('{} skip site low ac {} {}/{} {}'.format(fn, ac, c, n, lst[0:7]))
        continue
      # HWE
      data = gt2cnts(gts)
      hw, ea = HWtest(data[0], data[1], data[2])
      if ea < data[1] and hw < hwth:
        print('{} skip site HW ExHet {} {}|{} {}/{} {}'.format(fn, hw, data[1], round(ea, 3), c, n, lst[0:7]))
        continue
      for pop in ('nw6k', 'hd'):
        popi = [i for i, s in enumerate(samples) if spop[s] == pop]
        data = gt2cnts([gts[i] for i in popi])
        hw, ea = HWtest(data[0], data[1], data[2])
        if hw < hwth:
          print('{} skip site HW {} {} {}|{} {}/{} {}'.format(fn, pop, hw, data[1], round(ea, 3), c, n, lst[0:7]))
          break
      if hw < hwth: continue
      mhw, mpop, md = 0, '', []
      for pop in ('AFR', 'EAS', 'EUR', 'AMR', 'SAS'):
        popi = [i for i, s in enumerate(samples) if spop[s] == pop] # [ds[s] for s in samples if s in ds and spop[s] == pop]
        data = gt2cnts([gts[i] for i in popi])
        hw, ea = HWtest(data[0], data[1], data[2])
        if mhw < hw: mhw, mpop, md, mea = hw, pop, data, ea
      if mhw < hwth2:
        print('{} skip site HW 1KGP {} {} {}|{} {}/{} {}'.format(fn, mpop, mhw, md[1], round(mea, 3), c, n, lst[0:7]))
      c += 1
      yield ost + gts
  print('{} end of vcf {}/{}'.format(fn, c, n))

vis = [vcf_iter(fn) for fn in vcffiles]
vs = [next(vi) for vi in vis]

l = len(vcffiles)
endi = [0] * l

header = ''
if args.header is not None:
  for line in open(args.header):
    header += line

outfile = open(args.output, 'w')
for lst in io.splitIter(vcffiles[0]):
  if lst[0].startswith('##'): outfile.write(io.tabjoin(lst, '\n'))
  elif lst[0].startswith('#'):
    outfile.write(header)
    outfile.write(io.tabjoin(lst[0:9], samples, '\n'))
    break

while sum(endi) < l:
  ps = [v[1] for v in vs if v is not None]
  mp = min(ps)
  out = []
  for i in range(l):
    if endi[i] == 1: continue
    if vs[i][1] == mp:
      out.append(vs[i])
      try: vs[i] = next(vis[i])
      except StopIteration:
        vs[i] = None
        endi[i] = 1
  if len(out) > 0: out.sort()
  for v in out: outfile.write(io.tabjoin(v, '\n')) ## if no duplicates
  sys.stdout.flush()
       
