import gzip
import sys
import time
import argparse
from zbio import io 
#from scipy.stats import chi2
import numpy as np

parser = argparse.ArgumentParser()
#parser.add_argument('-i', dest='input', required=True, type=str, help='input vcf files with chr replaced by {}')
parser.add_argument('-o', dest='output', required=True, type=str, help='output stats file')
#parser.add_argument('-c', dest='chr', required=True, type=str, help='chr')
parser.add_argument('-t', dest='type', required=True, type=str, help='Variant type')
parser.add_argument('-S', dest='source', type=str, help='Variant source')
parser.add_argument('-s', dest='samples', type=str, default='rand100.txt', help='Selected samples')
#parser.add_argument('-p', dest='pops', type=str, default='../allpops.txt', help='Sample populations')
#parser.add_argument('--header', type=str, help='vcf headers to add')
parser.add_argument('--sparse', action="store_true", help='Use sparse AF intervals')
parser.add_argument('--sparse2', action="store_true", help='Use sparse2 AF intervals')
parser.add_argument('--keeptag', action="store_true", help='Keep tag variants')
parser.add_argument('--log', dest='log', type=str, help='output log dir')

args = parser.parse_args()

chrs = list(map(str, range(1, 23)))
afths = [0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
if args.sparse:
  afths = [0.005, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
elif args.sparse2:
  afths = [0.005, 0.025, 0.05, 0.2, 0.4, 0.6, 1.0]

an = 967 * 2.0
print('Job begin', time.ctime())
samples = []
for lst in io.splitIter(args.samples):
  samples.append(lst[0])

spop = {}
#for lst in io.splitIter(args.pops):
#  spop[lst[0]] = lst[2]
infos, skip = {}, {}
for chr in chrs:
  #for lst in io.splitIter('20240919zuihouyicimendal/high_error_rate_sites{}.txt'.format(chr)): ##
  #  key = chr, int(lst[1]), lst[2], lst[3]
  #  skip[key] = 1
  for lst in io.splitIter('../panel2/svtypes/svtypes_chr{}.txt.gz'.format(chr), skip=1):
  #for lst in io.splitIter('panel/chr{}.unrel.siteonly.vcf.gz'.format(chr)):
    #if lst[0].startswith('#'): continue
    vtype, source = lst[0], lst[1] # ''
    #if lst[2].startswith(chr):
    #  if len(lst[3]) == 1 and len(lst[4]) == 1: vtype = 'SNP'
    #  else: vtype = 'INDEL'
    #elif lst[2].startswith('phase2_'): vtype = 'ME'
    #elif lst[2] == '.':
    #  if lst[4].startswith('<'): vtype = 'SV'
    #  else: vtype = 'STR'
    #else: vtype = 'Unknown'
    #if vtype == 'SV' and args.type in ('INS', 'DEL', 'DUP', 'INV'):
    #  vtype = lst[4][1:-1]
    #elif vtype == 'ME' and args.type in ('MEI', 'MEA'):
    #  if lst[2].startswith('phase2_abs'): vtype = 'MEA'
    #  elif lst[2].startswith('phase2_ins'): vtype = 'MEI'
    if vtype != args.type: continue 
    if args.source is not None and source != args.source: continue
    
    key = lst[4], int(lst[5]), lst[6], lst[7], lst[8]
    if key in skip: continue
    ac, hetf = lst[2], float(lst[3]) # = an = '.'
    #for s in lst[7].split(';'):
    #  lst2 = s.split('=')
    #  if lst2[0] == 'AC': ac = lst2[1]
    #  elif lst2[0] == 'AN': an = lst2[1]
    #if lst[3] == '.': lst[3] = 'N'
    infos[key] = [vtype, source, ac, hetf]
    if len(infos) == 1: print(key, infos[key])
  print(chr, 'infos', len(infos), time.ctime())
#keptinfo = {'SVTYPE':1, 'SVLEN':1, 'END':1, 'MEI':1, '0START':1, '0END':1, 'MESTRAND':1}
#for fn in ('0830chuan/mea.vcf.gz','0830chuan/mei.vcf.gz', '0830chuan/SVchr{}.vcf.gz'.format(args.chr)):

keep = {}
#for chr in chrs:
#  for lst in io.splitIter('out/ldfilt.chr{}.txt'.format(chr)):
#    key = chr, int(lst[1]), lst[3], lst[4]
#    keep[key] = 1

gtc = {'0|0': 0, '0/0': 0, '0/1': 1, '1/0':1, '0|1': 1, '1|0': 1, '1|1': 2, '1/1': 2}
def vcf_iter(fn):
  #mac, miss, gqth, hwth, hwth2 = 2, 0.9, 5, 1e-6, 1e-10
  n = c = 0
  for lst in io.splitIter(fn):
    if lst[0].startswith('##'): continue
    if lst[0].startswith('#'):
      ds = {}
      for i in range(9, len(lst)): ds[lst[i]] = i
      si = [ds[s] for s in samples]
      #if None in si:
      #  for s in samples:
      #    if s not in ds: sys.stderr.write('{} sample not in vcf {}\n'.format(fn, s))
    else:
      ost = lst[0:8] # + ['GT']
      ost[0] = ost[0].replace('chr', '')
      #if ost[0] != args.chr: continue
      n += 1
      ost[1] = int(ost[1])
      if ost[3] == '.': ost[3] = 'N'
      key = ost[0], ost[1], ost[2], ost[3], ost[4]
      #if key not in keep: continue
      if key not in infos:
        #sys.stderr.write('Key not in infos:', key)
        continue
      #if infos[key][0] == 'STR': continue ####
      if ',' in ost[4]: continue #### multiallele
      #if key in infos: ost[7] = '{};{}'.format(infos[key], ost[7])
      gts0 = [lst[i].split(':')[0] for i in si]
      gts = [gtc[s] if s in gtc else None for s in gts0]
      fmt = lst[8].split(':')
      dss = None
      if 'DS' in fmt:
        qi = fmt.index('DS')
        dss = [float(lst[i].split(':')[qi]) for i in si]
      #l = len(si)
      #nmis = ac = 0
      c += 1
      yield ost + [gts, dss]
  sys.stderr.write('{} end of vcf {}/{}, {} {} {}\n'.format(fn, c, n, args.type, args.source, time.ctime()))

def dsr2(gts, dss):
  corr = np.corrcoef(gts, dss)[0][1]
  return corr * corr

def gtcon(gts1, gts2):
  l, c = len(gts1), 0.0
  for i in range(l):
    if gts1[i] == gts2[i]: c += 1
    elif abs(gts1[i] - gts2[i]) == 1: c += 0.5
  return c / l

#vis = [vcf_iter(fn) for fn in vcffiles]
#vs = [next(vi) for vi in vis]
#infn = 'sgdp/chr{}.full.vcf.gz'.format(args.chr)
#vi = vcf_iter(infn)
#v = next(vi)
infd = {'SNP': 'sgdp2/chr{}.snps.vcf.gz', 'INDEL': 'sgdp2/chr{}.indels.vcf.gz', 'ME': 'sgdp2/chr{}.svs.vcf.gz', 
        'MEI': 'sgdp2/chr{}.svs.vcf.gz', 'MEA': 'sgdp2/chr{}.svs.vcf.gz', 'SV': 'sgdp2/chr{}.svs.vcf.gz',
        'DUP': 'sgdp2/chr{}.svs.vcf.gz', 'INS': 'sgdp2/chr{}.svs.vcf.gz', 'DEL': 'sgdp2/chr{}.svs.vcf.gz',
        'INV': 'sgdp2/chr{}.svs.vcf.gz', 'STR': 'sgdp2/chr{}.str.vcf.gz',
}
#l = len(vcffiles)
#endi = [0] * l
l = len(samples)

#outfile = open(args.output, 'w')
#outfile = sys.stdout
#outfile.write('CHROM\tPOS\tID\tREF\tALT\tVType\tAC\tAN\tAC2\tAN2\tACI\tImputed\tR2DOS\tR2GT\tCC\n')
data = []
tags = {}
for s in samples:
  print(s, time.ctime())
  d = []
  for f in afths: d.append([[], [], [], []])  # true GT, DS, imp GT, AF
  data.append(d)
  for chr in chrs:
    infn = 'data/{0}/{0}.chr{1}.op.b5.vcf.gz'.format(s, chr) # infd[args.type].format(chr)
    for lst in io.splitIter(infn):
      if lst[0].startswith('#'): continue
      if ',' in lst[4]: continue
      key = lst[0], int(lst[1]), lst[2], lst[3], lst[4]
      if key not in infos: continue
      if key in tags: continue
      if not lst[7].endswith(';IMP'):
        tags[key] = 1
        continue
      ac = infos[key][2]
      af = float(ac) / an
      for j, ft in enumerate(afths):
        if af <= ft: break
      gt = gtc[lst[9].split(':')[0]] #
      fmt = lst[8].split(':')
      qi = fmt.index('DS')
      ds = float(lst[9].split(':')[qi])
      d[j][1].append(ds)
      d[j][2].append(gt)
      #d[j][3].append(af)

for chr in chrs:
  #vi = vcf_iter(infn)
  #v = next(vi)
  for p in vcf_iter('../panel2/chr{}.filt.vcf.gz'.format(chr)): # args.input.format(chr)):
    #while v is not None and v[1] <= p[1]:
    #  if (v[1], v[3], v[4]) == (p[1], p[3], p[4]): break
    #  try: v = next(vi)
    #  except StopIteration:
    #    v = None
    #    break
    #if v is None or v[1] > p[1]:
      #print('Var not found in full data: {}'.format(p[0:5]))
    #  continue
  #print(len(v[8]), len(p[9]), v[8], p[9])
    #typed = 'T'
    #if p[7].startswith('IMPUTED'): typed = 'I'
    #elif p[7].endswith(';IMP'): typed = 'I'
    #if not args.keeptag and typed == 'T': continue

    key = p[0], p[1], p[2], p[3], p[4]
    if key in tags: continue
    ac = infos[key][2]
    af = float(ac) / an
    for j, ft in enumerate(afths):
      if af <= ft: break
    for i in range(l):
      #if v[8][i] is None: continue
      d = data[i][j]
      d[0].append(p[8][i])
      #d[1].append(p[9][i])
      #d[2].append(p[8][i])
      d[3].append(af)
    #r2 = dsr2(v[8], p[9])
    #r2g = dsr2(v[8], p[8])
    #con = gtcon(v[8], p[8])
    #key = p[1], p[2], p[3], p[4]
    #if key not in infos: vinfo = ['.', 0, an]
    #else:
    #  vinfo = infos[key]
    #  an = vinfo[2]
    #ac2, an2, aci = sum(v[8]), len(v[8]) * 2, sum(p[8])
    #for s in p[7].split(';'):
    #  lst2 = s.split('=')
    #  if lst2[0] in ('DR2', 'R2'):
    #    aci = '{}|{}'.format(aci, lst2[1])
    #    break
    #outfile.write(io.tabjoin(v[0:5], vinfo, ac2, an2, aci, typed, r2, r2g, con, '\n'))

outfile = open(args.output, 'w')
header = ['Cnt_{0}\tAF_{0}\tR2DS_{0}\tR2GT_{0}\tCC_{0}'.format(str(f * 100).replace('.0', '')) for f in afths]
outfile.write(io.tabjoin('Sample', header, '\n'))
for i in range(l):
  out = [samples[i]]
  for j, ft in enumerate(afths):
    d = data[i][j]
    n = len(d[0])
    if n == 0: af = r2d = r2g = cc = 'NA'
    else:
      af = round(sum(d[3]) / n, 5)
      #print(i, j, af) #, d[0], d[1])
      r2d = round(dsr2(d[0], d[1]), 5)
      r2g = round(dsr2(d[0], d[2]), 5)
      cc = round(gtcon(d[0], d[2]), 5)
    out += [n, af, r2d, r2g, cc]
  outfile.write(io.tabjoin(out, '\n'))

print('Done', time.ctime())


