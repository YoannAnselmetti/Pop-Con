###
###   Author:              Yoann Anselmetti
###   Last modification:   2019/08/06
###
###   License: This software is distributed under the CeCILL free software license (Version 2.1 dated 2013-06-21)
###

from os import path, makedirs
import errno


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc:
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise



