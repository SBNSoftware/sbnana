class SelectionInfo:

  def __init__(self, Name, Cut, SpillCut="kNoSpillCut", Comment=""):
    self.Name = Name
    self.Cut = Cut
    self.SpillCut = SpillCut
    self.Comment = Comment

  def GetSampleSelectionLines(self):

    return '''  //==== %s
  baseSampleNames.push_back("%s");
  baseSampleCuts.push_back(%s);
  baseSampleSpillCuts.push_back(%s);
'''%(self.Comment, self.Name, self.Cut, self.SpillCut)

  def GetEventSelectionLines(self):

    return '''  //==== %s
  cutNames.push_back("%s");
  cuts.push_back(%s);
  spillcuts.push_back(%s);
'''%(self.Comment, self.Name, self.Cut, self.SpillCut)
