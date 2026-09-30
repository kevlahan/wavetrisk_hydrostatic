"""Check native soil arithmetic against the legacy kernels after storage-only normalization."""
import hashlib,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
NAMES=['Compute_scalar_wavelets','Compute_velo_wavelets','Compute_velo_wavelets_penta','velo_interp_penta_corr',
       'Interp_inner_velo','Interp_outer_velo','Restrict_velo','Interp_node','compress_scalar','compress_vector',
       'Prolong_scalar','Reconstruct_scalar','Reconstruct_outer_velo','Reconstruct_inner_velo','Reconstruct_velo_penta']
def routine(source,name):
    m=re.search(r'^  (?:real\(dp\) )?(?:subroutine|function) '+name+r'\b.*?^  end (?:subroutine|function) '+name+r'\b',source,re.M|re.S|re.I)
    assert m,name
    return m[0]
def normalize(source,native=False):
    s=' '.join(line.split('!')[0] for line in source.splitlines()).replace('&','')
    if native:
        while m:=re.search(r'call\s+set_(soil_\w+)\s*\(',s,re.I):
            start=m.end();j=start;depth=1;comma=None
            while depth:
                if s[j]=='(':depth+=1
                elif s[j]==')':depth-=1
                elif s[j]==',' and depth==1 and comma is None:comma=j
                j+=1
            assert comma is not None
            s=s[:m.start()]+m[1]+'('+s[start:comma]+')='+s[comma+1:j-1]+s[j:]
        for name in NAMES:s=re.sub(r'\bsoil_'+name+r'\b',name,s,flags=re.I)
        for a,b in [('soil_scalar_wavelet','wc_s'),('soil_vector_wavelet','wc_u'),('soil_scalar','scalar'),('soil_vector','velo'),('base_weight','Iu_Base_Wgt')]:
            s=re.sub(r'\b'+a+r'\b',b,s,flags=re.I)
    return re.sub(r'\s+','',s).lower()
wave=(BASE/'src/wavelet.f90').read_text();adapt=(BASE/'src/adapt.f90').read_text();native=(BASE/'src/parallel_block_inverse.f90').read_text();rows={}
def base_weight(source,name):
    text=source.replace('&','')
    m=re.search(r'\b'+name+r'\(9\)\s*=\s*(\[.*?\]\s*/\s*16\.0_dp)',text,re.S|re.I)
    assert m,name
    return re.sub(r'\s+','',m[1]).lower()
assert base_weight(wave,'Iu_Base_Wgt')==base_weight(native,'base_weight')
for name in NAMES:
    old=normalize(routine(adapt if name.startswith('compress_') else wave,name))
    new=normalize(routine(native,'soil_'+name),True)
    assert old==new,name
    rows[name]=hashlib.sha256(old.encode()).hexdigest()
save(BASE/'docs/parallel-blocks-stage247-kernels.json',dict(status='passed',kernels=rows,
     normalization='Ignore comments/whitespace/continuations; rename native storage accessors and expand setter calls back to assignments.',
     sources={name:digest(BASE/name) for name in ['src/wavelet.f90','src/adapt.f90','src/parallel_block_inverse.f90']}))
print('PASS fifteen native soil arithmetic kernels',flush=True)
