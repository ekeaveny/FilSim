# Copies seed par, bak and springlinks.dat files from one folder to another.
# Takes two arguments.

from __future__ import division, print_function
import os
import sys
import shutil

if len(sys.argv) < 3:
    print("ERROR: You need to provide two arguments: the folder (timestamp) "
          "you are copying from, and the folder you are copying to.")
    sys.exit()

folder_from = sys.argv[1]
folder_to = sys.argv[2]

start_dir = os.getcwd()

os.chdir(folder_from)

subdirectories_from = [name for name in os.listdir(".") if os.path.isdir(name)]

print ("Searching source folder...")

files_from = [[] for s in subdirectories_from]
for s, subdirectory in enumerate(subdirectories_from):
    os.chdir(subdirectory + "/output")
    data_files = [name for name in os.listdir(".") if os.path.isfile(name)]
    os.chdir('../..')

    backup_files = [name for name in data_files if name[-4:] == '.bak']
    parameter_files = [name for name in data_files if name[-4:] == '.par']
    springlinks_files = [name for name in data_files if name[-16:] == '-springlinks.dat']

    is_a_spring_sim = len(springlinks_files) == 1
    files_are_there_that_i_expect = (len(backup_files) == 1
                                     and len(parameter_files) == 1)

    if files_are_there_that_i_expect:
        files_from[s].append(backup_files[0])
        files_from[s].append(parameter_files[0])

        if is_a_spring_sim:
            files_from[s].append(springlinks_files[0])

os.chdir(start_dir)
os.chdir(folder_to)

subdirectories_to = [name for name in os.listdir(".") if os.path.isdir(name)]

confirmed_files_from = []
confirmed_files_to = []

print ("...done.")
print ()
print ("Searching destination folder...")

overwrite_warnings = 0
for s, subdirectory in enumerate(subdirectories_from):
    if subdirectory in subdirectories_to:
        # Is there a file of this name already in the dest folder?
        os.chdir(subdirectory + "/output")
        data_files = [name for name in os.listdir(".") if os.path.isfile(name)]
        os.chdir('../..')
        for file_old in files_from[s]:
            if file_old in data_files:
                if file_old[-4:] != ".par":
                    print("(!) " + file_old + " already exists in "
                          + folder_to + ". File will not be overwritten.")
                    overwrite_warnings += 1
                    files_from[s].remove(file_old)

        from_directory = folder_from + "/" + subdirectory + "/output/"
        files_from_full_path = [from_directory + i for i in files_from[s]]
        files_to_full_path = folder_to + "/" + subdirectory + "/output/"

        confirmed_files_from = confirmed_files_from + files_from_full_path
        confirmed_files_to = (confirmed_files_to
                              + [files_to_full_path for i in files_from_full_path])

print ("...done.")
print ()
print("Copying...")

os.chdir(start_dir)
f = -1
for f, src in enumerate(confirmed_files_from):
    dst = confirmed_files_to[f]
    print(str(f + 1) + ". " + src + "\n--> " + dst)
    print()
    shutil.copy2(start_dir + "/" + src, start_dir + "/" + dst)

print("...done.")

print("Finished. " + str(f + 1) + " files copied.")
if overwrite_warnings > 0:
    print(str(overwrite_warnings) + " files were not copied because they would "
          "overwrite existing files (see log above).")
