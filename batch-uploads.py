# Code for uploading

import os
directory = 'output/2003091044a'

for filename in os.listdir(directory):
	name_no_ext = filename[:-4]
	print(filename)
	os.system('scp "%s" "%s:%s"' % (directory+"/"+filename, "cx2", "~/MultiFilamentFCM-Adam/2003261341/"+name_no_ext+"/output/") )

# Code for downloading
# scp 'rizzuto:~/MultiFilamentFCM-Adam/2003091044/*/output/*.bak' ./2003091044a/
