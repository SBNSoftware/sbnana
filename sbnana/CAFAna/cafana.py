if __name__ == '__main__':
    print("This is the source file that implements the cafana python module.")
    print("You should 'import cafana' at the top of a .py script")
    print("and execute it with 'cafe myscript.py'")
    exit(1)

import os

if 'SBNANA_INC' in os.environ:
    inc = os.environ['SBNANA_INC']
else:
    inc = os.environ['MRB_INSTALL']+'/sbnana/'+os.environ['SBNANA_VERSION']+'/include/'

if 'SBNANAOBJ_INC' in os.environ:
    sao_inc = os.environ['SBNANAOBJ_INC']
else:
    sao_inc = os.environ['MRB_INSTALL']+'/sbnanaobj/'+os.environ['SBNANAOBJ_VERSION']+'/include/'

os.environ['ROOT_INCLUDE_PATH'] = \
  ':'.join([inc,
            inc+'sbnana',
            inc+'sbnana/CAFAna',
            sao_inc,
            os.environ['EIGEN_INC'],
            os.environ['OSCLIB_INC'],
            os.environ['SRPROXY_INC']])

import ROOT

if 'SBNANA_FQ_DIR' in os.environ:
    ROOT.gApplication.ExecuteFile('${SBNANA_FQ_DIR}/bin/rootlogon.C')
else:
    ROOT.gApplication.ExecuteFile('$MRB_BUILDDIR/sbnana/bin/rootlogon.C')
print('  in python')
ROOT.gROOT.ForceStyle()

print('Load libraries...')
for lib in ['Minuit2',
            'sbnanaobj_StandardRecordProxy',
            'CAFAnaCore',
            'SBNAnaVars',
            'SBNAnaCuts',
            'CAFAnaSysts',
            'CAFAnaExtrap',
            'CAFAnaPrediction',
            'CAFAnaExperiment',
            'CAFAnaAnalysis',
            'SBNAnaVars',
            'SBNAnaCuts']:
    print(' ', lib)
    ROOT.gSystem.Load('lib'+lib+'.so')


import cppyy

print('Load dictionaries... (please ignore errors about .pcm files)')
for d in ['CAFAna', 'SBNAna']:
    print(' ', d)
    cppyy.load_reflection_info('lib'+d+'_dict.so')

class PyCAFAna:
    def __init__(self, cppyy):
        self._cppyy = cppyy

    def CSliceVar(self, body):
        '''Construct a new slice Var given the C++ body as a string'''
        var = 'pyvar_'+self._cppyy.gbl.ana.UniqueName()
        text = '#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"\ndouble '+var+'_func(const caf::SRSliceProxy* srp){\nconst caf::SRSliceProxy& sr = *srp;\n'+body+'\n}\nconst ana::Var '+var+'('+var+'_func);'
        self._cppyy.cppdef(text)
        return getattr(self._cppyy.gbl, var)

    def CSpillVar(self, body):
        '''Construct a new spill Var given the C++ body as a string'''
        var = 'pyvar_'+self._cppyy.gbl.ana.UniqueName()
        text = '#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"\ndouble '+var+'_func(const caf::SRSpillProxy* srp){\nconst caf::SRSpillProxy& sr = *srp;\n'+body+'\n}\nconst ana::Var '+var+'('+var+'_func);'
        self._cppyy.cppdef(text)
        return getattr(self._cppyy.gbl, var)

    def SimpleSliceVar(self, name):
        '''Equivalent of the SIMPLEVAR() macro'''
        return self.CSliceVar('return sr.'+name+';')

    def SimpleSpillVar(self, name):
        '''Equivalent of the SIMPLEVAR() macro'''
        return self.CSpillVar('return sr.'+name+';')

    def CSliceCut(self, body):
        '''Construct a new Cut given the C++ body as a string'''
        cut = 'pycut_'+self._cppyy.gbl.ana.UniqueName()
        text = '#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"\nbool '+cut+'_func(const caf::SRSliceProxy* srp){\nconst caf::SRSliceProxy& sr = *srp;\n'+body+'\n}\nconst ana::Cut '+cut+'('+cut+'_func);'
        self._cppyy.cppdef(text)
        return getattr(self._cppyy.gbl, cut)

    def CSpillCut(self, body):
        '''Construct a new Cut given the C++ body as a string'''
        cut = 'pycut_'+self._cppyy.gbl.ana.UniqueName()
        text = '#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"\nbool '+cut+'_func(const caf::SRSpillProxy* srp){\nconst caf::SRSpillProxy& sr = *srp;\n'+body+'\n}\nconst ana::Cut '+cut+'('+cut+'_func);'
        self._cppyy.cppdef(text)
        return getattr(self._cppyy.gbl, cut)


    # If we don't provide it explicitly, assume it's from the ana namespace
    def __getattr__(self, name):
        return getattr(self._cppyy.gbl.ana, name)

import sys
sys.modules['cafana'] = PyCAFAna(cppyy)
