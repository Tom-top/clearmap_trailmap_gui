# -*- coding: utf-8 -*-
"""
Module providing some sound output to signal processes are done

Example:
  >>> import ClearMap.Utils.Sound as snd
  >>> snd.beep()
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'

import os


def beep(duration = 1, frequency = 340):
  os.system('(speaker-test -t sine -f %d >/dev/null)& pid=$! ; sleep 1s ; kill -9 $pid' % (frequency))
  os.system('echo -e "\a" >/dev/null');
  
 
if __name__ == "__main__":
  import ClearMap.Utils.Sound as snd;
  reload(snd)
  snd.beep(frequency= 300)  