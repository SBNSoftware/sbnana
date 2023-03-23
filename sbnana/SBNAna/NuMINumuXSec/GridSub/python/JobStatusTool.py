import os
import Logger

class JobStatusTool:

  def __init__(self, UserName, RanStr, LogThreshold="INFO"):

    self.logger = Logger.Logger("JobStatusTool", LogThreshold)
    self.UserName = UserName
    self.RanStr = RanStr

  def GetJobStatusByUser(self):

    filepath_JobStatus = ('/tmp/jobstatus_%s_%s.log')%(self.RanStr, self.UserName)

    os.system('jobsub_q -G icarus --user %s &> %s'%(self.UserName, filepath_JobStatus))

    failed = False
    for l in open(filepath_JobStatus).readlines():
      if "Failed to fetch ads" in l:
        failed = True
        break

    if failed:
      self.GetJobStatusByUser()

    return ('/tmp/jobstatus_%s_%s.log')%(self.RanStr, self.UserName)

  def GetJobStatusByJobID(self, lines_JobStatus, this_jobID):
    for line in lines_JobStatus:

      # '33667526.0@jobsub03.fnal.gov          jskim           06/13 19:53   0+00:30:20 H   0   0.0 test.sh_20220613_195324_2194810_0_1_wrap.sh \n'
      words = line.split()
      '''
      0: "33667526.0@jobsub03.fnal.gov"
      1: "jskim"
      2: "06/13"
      3: "19:53"
      4: "0+00:30:20"
      5: "H"
      '''

      self.logger.LOG("DEBUG", words, "GetJobStatusByJobID")

      if len(words)==0:
        break

      if words[0]==this_jobID:
        jobFlag = words[5]
        if jobFlag=="H":
          return "HOLD"
        elif jobFlag=="R":
          return "RUNNING::%s"%(words[4])
        elif jobFlag=="I":
          return "IDLE"
        elif jobFlag=="C":
          return "FINISHED"
        else:
          return jobFlag
       
    return "FINISHED"
