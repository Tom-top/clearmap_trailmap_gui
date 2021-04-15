# -*- coding: utf-8 -*-
"""
Utility module to convert a regular expression to patterns and back

See the 'sre_parse' python library for definition of the pattern
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__docformat__ = 'rest'


# invert regular expression

import sre_parse as sre

TEXT = 'text';

class PatternInverter():
  # inverted categories
  inCategories =    { 
      sre.CATEGORY_DIGIT :  r"\d", 
      sre.CATEGORY_NOT_DIGIT: r"\D", 
      sre.CATEGORY_SPACE : r"\s",
      sre.CATEGORY_NOT_SPACE : r"\S",
      sre.CATEGORY_WORD : r"\w",
      sre.CATEGORY_NOT_WORD : r"\W"
      };
  
  atCategories = {
      sre.AT_BEGINNING_STRING : r"\A", 
      sre.AT_BOUNDARY : r"\b",
      sre.AT_NON_BOUNDARY : r"\B", 
      sre.AT_END_STRING : r"\Z",
      sre.AT_BEGINNING : "^",
      sre.AT_END : "$"
      };
  
  escapes = {
      ord("\a"): r"\a",
      ord("\b"): r"\b",
      ord("\f"): r"\f",
      ord("\n"): r"\n",
      ord("\r"): r"\r",
      ord("\t"): r"\t",
      ord("\v"): r"\v",
      ord("\\"): r"\\"
      }
  
  def __init__(self, groups = None):
    self.groups = groups;
  
  def generateAny(self, pattern):
    return '.';
  
  def generateLiteral(self, pattern):
    _, val = pattern; 
    c = self.escapes.get(val);
    if c is not None:
      return c;
    c = chr(val);
    if c in sre.SPECIAL_CHARS:
      c = '\\' + c;
    return c
    
  def generateNotLiteral(self, pattern):
    _, val = pattern; 
    c = self.escapes.get(val);
    if c is not None:
      return "^" + c;
    c = chr(val);
    if c in sre.SPECIAL_CHARS:
      c = '\\' + c;
    return "^" + c;
  
  def generateMaxRepeat(self, pattern):
    rmin = pattern[1][0];
    rmax = pattern[1][1];
    rep = pattern[1][2];
    
    if rmin == 0 and rmax == sre.MAXREPEAT:
      return self.generateRE(rep) + '*';
    if rmin == rmax:
      return self.generateRE(rep) + '{' + str(rmin) + '}';
    return self.generateRE(rep) + '{' + str(rmin) + ',' + str(rmax) + '}';
    
  def generateMinRepeat(self, pattern):
    rep = pattern[1][2];
    return self.generateRE(rep) + '*?';
    
  def generateGroup(self, pattern, groupname):
    return '(?P<' + groupname + '>' + self.generateRE(pattern) + ')';
    
  def generateGroupRef(self, pattern):
    gid = pattern[1];
    if gid in self.groups.keys():
      return '(P?=' + self.groups[gid] + ')';
    else:
      return '\\' + str(gid)
   
  def generateNegate(self, pattern):
    return '^';
  
  def generateRange(self, pattern):
    r = pattern[1];
    return chr(r[0]) + '-' + chr(r[1]);
  
  def generateIn(self, pattern):
    p = pattern[1];
    if len(p) == 1:
      tp = p[0];
      if tp[0] == sre.CATEGORY:
        return self.inCategories[tp[1]];
      elif tp[0] == sre.AT:
        return self.atCategories[tp[1]];
    return '[' + self.generateRE(p) + ']';
  
  def generateBranch(self, pattern):
    return '|'.join([self.generateRE(p) for p in pattern[1][1]]);
  
  def generateSubpattern(self, pattern):    
    sid = pattern[1][0];
    if self.groups is not None and sid in self.groups.keys():
      groupname = self.groups[sid];
      return self.generateGroup(pattern[1][1], groupname);
    
    if sid is not None:
      return '(' + self.generateRE(pattern[1][1]) + ')';
    
    sp = pattern[1][1];
    if sp[0][0] == sre.GROUPREF_EXISTS:
      gr = sp[0];
      gid = gr[1][0];
      gyes = gr[1][1];
      gno = gr[1][2];
      if self.groups is not None and gid in self.groups.keys():
        groupname = self.groups[gid];
      else:
        groupname = str(gid);
      if gno is None:
        return '(?(' + groupname + ')' + self.generateRE(gyes) + ')';
      else:
        return '(?(' + groupname + ')' + self.generateRE(gyes) + '|' + self.generateRE(gno) + ')';
    else:
      return '(?:' + self.generateRE(pattern[1][1]) + ')';
  
  def generateAssert(self, pattern):    
    return '(?=' + self.generateRE(pattern[1][1]) + ')';
    
  def generateAssertNot(self, pattern):    
    return '(?!' + self.generateRE(pattern[1][1]) + ')';
   
  def generateAt(self, pattern):
    return self.atCategories[pattern[1]];
  
  def generateREType(self, pattern):
    ptype = pattern[0];
    if ptype == sre.LITERAL:
      return self.generateLiteral(pattern);
    if ptype == sre.NOT_LITERAL:
      return self.generateNotLiteral(pattern);
    elif ptype == sre.ANY:
      return self.generateAny(pattern);
    elif ptype == sre.AT:
      return self.generateAt(pattern);
    elif ptype == sre.ASSERT:
      return self.generateAssert(pattern);
    elif ptype == sre.ASSERT_NOT:
      return self.generateAssertNot(pattern);
    elif ptype == sre.SUBPATTERN:
      return self.generateSubpattern(pattern);
    elif ptype == sre.BRANCH:
      return self.generateBranch(pattern);
    elif ptype == sre.IN:
      return self.generateIn(pattern);   
    elif ptype == sre.GROUPREF:
      return self.generateGroupRef(pattern);  
    elif ptype == sre.MIN_REPEAT:
      return self.generateMinRepeat(pattern);
    elif ptype == sre.MAX_REPEAT:
      return self.generateMaxRepeat(pattern);
    elif ptype == sre.RANGE:
      return self.generateRange(pattern);
    elif ptype == sre.NEGATE:
      return self.generateNegate(pattern); 
    elif ptype == TEXT:
      return pattern[1];
    else:
      raise RuntimeError("Inversion of pattern of type '%s' not implemented!" % ptype);    
  
  def generateRE(self, pattern):
    r = '';    
    for p in pattern:
      #print p
      r += self.generateREType(p);
    return r;    
 

def expressionToPattern(expression):
  """Convert regular expression to a parsed pattern for manipulation
  
  Arguments:
    expression (str): regular expression
    
  Returns:
    object: parsed pattern
  """
  return sre.parse(expression);


def patternToExpression(pattern):
  """Convert pattern to regular expression
  
  Arguments:
    pattern (object): regular expression pattern
    
  Returns:
    str: regular expression
  """  
  
  if isinstance(pattern, sre.SubPattern) or isinstance(pattern, sre.Pattern):
    gd = pattern.pattern.groupdict;
  else:
    gd = {};
  
  gdi = {y:x for x,y in gd.items()};
  pi = PatternInverter(groups = gdi);
  return pi.generateRE(pattern);


def test():
  import ClearMap.Utils.InverseRegularExpression as ire;
  import sre_parse as sre;

  reload(ire)  
  source = '/test/test_(?P<row>\d{4})_(?P<col>\d{3}).tif';
  p = sre.parse(source);
  ire.patternToExpression(p)

  reload(ire)  
  source = r'/test/test_(?:\d)_(?P<col>\d{3})_[7-9][.](?=col)tif$';
  p = sre.parse(source);
  ire.patternToExpression(p)


if __name__ == "__main__":
  test();


