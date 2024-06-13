#!/opt/local/bin/python
# -*- coding: utf-8 -*-

########################
#
# this program run the whole program
#
########################

# 2018-12-29 @ Menaka
import main_code
try:
  main_code.main_act()
except Exception as e:
  print (e)


## check before run ################
#
# 1. compile CaMa-Flood
# 2. set input runoff data
# 3. set parameters in params.py