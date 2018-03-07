import os
import sys
import ROOT
from math import sqrt,pi

#-------------------------------------------------------------------------------
# functions

# write to file
def write_output():
    mixture_count = 0;
    # loop over tree, print variables
    for event in tree:
    #    if tree.size_nc_weight[1] == 12.0:
        if tree.ic_para == -1:
            # Event Number
            f.write(str(tree.size_nc_weight[0]) + "\n")
            # Compontent Number
            f.write(str(tree.size_nc_weight[1]) + "\n")
            # Tracking Particle 5er Vector
            for i in range(0,5):
                f.write(str(format(tree.tp_track[i], '.10f')) + "\n")
            # Indicater for total mixture = -1
            f.write(str(tree.ic_para) + "\n")
            # Weight of total mixture = 1
            f.write(str(format(tree.size_nc_weight[2], '.8f')) + "\n")
            # Assoc 5er Vector for total mixture
            for j in range(0,5):
                f.write(str(format(tree.localPars[j], '.10f')) + "\n")
            # Assoc 5x5 cov matrix for total mixture
            for a in range(0,5):
                for b in range(0,a+1):
                    f.write(str(format(tree.localCov(a,b), '.14f')) + "\n")
            mixture_count +=1
    #        print "\n"
            if (mixture_count % 500 == 0):
                print "mixture_count = " + str(mixture_count)
    #        print "total mixture written (" + str(tree.ic_para) + ")"
        if tree.ic_para > -1:
            # component index
            f.write(str(tree.ic_para) + "\n")
            # component weight
            f.write(str(format(tree.size_nc_weight[2], '.8f')) + "\n")
            # Assoc 5er Vector of current mixture
            for j in range(0,5):
                f.write(str(format(tree.localPars[j], '.10f')) + "\n")
            # Assoc 5x5 cov matrix of current mixture
            for a in range(0,5):
                for b in range(0,a+1):
                    f.write(str(format(tree.localCov(a,b), '.14f')) + "\n")
    #        print "mixture " + str(tree.ic_para) + " written"

    return mixture_count

#-------------------------------------------------------------------------------
# main code


# open root file & tree
pathname = "/afs/cern.ch/work/c/cbernkop/condor_output/output_ex1_180208_144044"
# remove file
try:
    os.remove(pathname + "/output.csv")
except:
    print("no output.csv file to remove found")

# loop over associator files
total_mixture_count = 0

i = 0
max_files = 1
while i <= max_files:
    # open csv, append to end of file
    f = open(pathname + "/output.csv","a")

    i += 1
# for opening files with associator executed by Christoph Bernkopf
   # filename = pathname + "/output_gsf_associator_" + str(i) + ".root"

# for opening files with associator executed by Wolfgang Adam
   # filename = "/afs/cern.ch/work/a/adamwo/Bernkopf/MultiElectron_FlatPt5To100_ext1/180126_081601/output.root"
    filename = "/afs/cern.ch/work/a/adamwo/Bernkopf/MultiElectron_FlatPt5To100_ext1/180208_144044/output.root"

# try to open new tree
    try:
        tf = ROOT.TFile(filename)
        tree_dir = tf.Get("TrackRecSimByDR")
#        tree_dir = tf.Get("MyTrackAssociator")
        tree = tree_dir.Get("track_associator_tree")

#        f = open(pathname + "/output_" + str(i) + ".csv","w+")

# write tree content into opened csv
        mixture_count = write_output()
        
        total_mixture_count += mixture_count
        print "tota_mixture_count = " + str(total_mixture_count)


# if error - stop while loop
    except:
        print(filename, "doesn't exist")
        continue

f.close()

print "Python finished."

f = open("output.csv","a")
f.truncate()
f.close()

