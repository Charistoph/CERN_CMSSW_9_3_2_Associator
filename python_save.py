import os
import sys
import ROOT
from math import sqrt,pi

#-------------------------------------------------------------------------------
# functions

# write to file
def write_output(path,iteration,filename):
    f = open(pathname + "/" + filename + "_" + iteration + ".csv","w+")

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

    f.close()

# create output folder
#outputpath = filename.replace(".root", "_csv")
##namestr = filename.replace(".root", "")
##outputpath = os.path.join('output_' + namestr)
#if (os.path.exists(outputpath)==False):
#    os.makedirs(outputpath)
#    print "\npath created", outputpath
#
#print "\nPath and file routine complete.\n"


#-------------------------------------------------------------------------------
# main code

# open root file & tree
pathname = "/afs/cern.ch/work/c/cbernkop/condor_output/assoc_output_20171120/"
os.chdir(pathname)

i = 0
while i < 100:
    i += 1
    path = pathname + str(i)
    filename = path + "/output_gsf_associator.root"

    try:
        tf = ROOT.TFile(filename)
        tree_dir = tf.Get("MyTrackAssociator")
        tree = tree_dir.Get("track_associator_tree")

        write_output(pathname,i,"output")

    except:
        break


print "Python finished."

