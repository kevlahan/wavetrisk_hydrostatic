"""Isolate campaign build settings without changing the user's shell environment."""
from contextlib import contextmanager
from pathlib import Path
import json,os
from common import save

# These inputs can alter make evaluation, flags, source selection or wrapper
# configuration. The sealed Makefiles/command line supply their intended values.
CONTROLLED_KEYS = (
 'FFLAGS','F90FLAGS','FCFLAGS','F77FLAGS','CPPFLAGS','CFLAGS','CXXFLAGS',
 'LDFLAGS','LDLIBS','LIBS','MAKEFLAGS','MFLAGS','GNUMAKEFLAGS','MAKEFILES',
 'MAKEOVERRIDES','MPIF90','FC','F90','F77','CC','CXX','LD','FFRONT',
 'DEBUG','FP_CONTRACT','STAGE202','STAGE203','TEST_CASE','PARAM','PREFIX',
 'SRC_DIR','BIN_DIR','BUILD_DIR','MODDIR','NF_CONFIG','NC_CONFIG',
 'LAPACK','NETLIB_LAPACK_ROOT','MPICH_FC','MPICH_F90','MPICH_F77',
 'MPICH_FCFLAGS','MPICH_F90FLAGS','MPICH_F77FLAGS','MPICH_CPPFLAGS','MPICH_LDFLAGS')

def normalized_environment(inherited,expected):
 clean=dict(inherited)
 removed={k:clean.pop(k) for k in CONTROLLED_KEYS if k in clean}
 pinned=dict(expected['build_dependencies'])
 if pinned['LAPACK']!='-llapack':raise ValueError('Unexpected sealed LAPACK linkage')
 clean.update(pinned)
 return clean,dict(inherited_controlled=removed,effective_pinned=pinned,
  policy='Discard inherited build overrides in this Python process only; sealed Makefiles, command line and job55551 dependency paths control compilation.')

@contextmanager
def controlled_environment(kit,out):
 original=dict(os.environ)
 expected=json.loads((Path(kit)/'expected-toolchain.json').read_text())
 clean,record=normalized_environment(original,expected)
 save(Path(out)/'build-environment313.json',record)
 try:
  os.environ.clear();os.environ.update(clean)
  yield record
 finally:
  os.environ.clear();os.environ.update(original)
