# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import urllib2
import codecs
import re
response = urllib2.urlopen("https://www.aqistudy.cn/historydata/monthdata.php?city=%E5%92%B8%E9%98%B3")
str = response.read()
pattern = re.compile('<td align="center">(.*?)</td>',re.S)
items = re.findall(pattern,str)
if(len(items) % 6 != 0):
    raise Exception("The length of items can not be divisible by 6 , there is an error !")
## 写入到文本文件pm25.txt
f2 = codecs.open("pm25.txt" , "w")
for i in range(0 , len(items) , 6):
    f2.write("gyk"+items[i+3] + "\n")
f2.close