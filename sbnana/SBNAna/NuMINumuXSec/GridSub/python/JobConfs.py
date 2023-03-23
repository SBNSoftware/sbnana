class JobConfs:

  def __init__(self, ana, ss, es, i, o, n):
    self.ana = ana
    self.ss = ss
    self.es = es
    self.i = i
    self.o = o
    self.n = n

  def GetCommand(self):

    cmd = 'create-batch'
    cmd += ' -a %s \\\n'%(self.ana)
    cmd += ' --ss %s \\\n'%(self.ss)
    cmd += ' --es %s \\\n'%(self.es)
    if 'txt' in self.i:
      cmd += ' -l %s \\\n'%(self.i)
    else:
      cmd += ' -i %s \\\n'%(self.i)
    cmd += ' -o %s \\\n'%(self.o)
    cmd += ' -n %d'%(self.n)

    return cmd
