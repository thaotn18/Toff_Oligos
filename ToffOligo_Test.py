from __future__ import print_function

import getopt, sys

import math
import warnings
import numpy


#from Bio.SeqUtils import MeltingTemp as mt
#from Bio.Seq import Seq
#from Bio import BiopythonWarning


NN_table_desc = ""
INN_table_desc = ""
TNN_table_desc = ""
DE_table_desc = ""



from Bio import SeqUtils, Seq
from Bio import BiopythonWarning
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

def GC(sequence):
    return 100 * gc_fraction(sequence, ambiguous="remove")

gc_fraction("ATT")


def disp_fork(ns,O5,Fork,O3):
    n = ns - Fork
    if O5 < Fork : O5 = Fork
    if O5 > ns-1 : O5 = ns-1
    if O3 < O5 : O3 = O5
    if O3 > ns : O3 = ns
    print("Fork = %d; O5 = %d; O3 = %d" %(Fork,O5,O3))
    for lin in range (0,ns+4):
        #line = ("%d %d " % (lin,n))
        line = (" ")
        if (lin == 0) :
            for col in range (0,5+(2*n)):
                line = line + " "
            if (O5 > 0):
                line = line + " "
                for col in range (0,O5-1):
                    line = line + "__"
                for col in range (O5-1,O5):
                    line = line + " "
            for col in range (2*O5,2*O3-1):
                line = line + " "
            if (O3 < n-1):
                line = line + " _"
                for col in range (O3+1,n-1):
                    line = line + "__"
                line = line + "_ "
        if (lin == 1) :
            for col in range (0,3+(2*n)):
                line = line + " "
            line = line + "5'"
            if (O5 > 0):
                for col in range (0,O5-1):
                    line = line + "/ "
                for col in range (O5-1,O5):
                    line = line + "/\\"
            for col in range (2*O5,2*O3-1):
                line = line + "_"
            if (O3 < ns):
                line = line + "/\\"
                for col in range (O3+1,ns):
                    line = line + " \\"
            line = line + "3'"
        if (lin == 2) :
            for col in range (0,5+(2*n)):
                line = line + " "
            for col in range (0,O5):
                line = line + "  "
            for col in range (O5,O3):
                line = line + "| "
            for col in range (O3,n+1):
                line = line + "  "
        if (lin == 3) :
            line = line + " "
            for col in range (n,ns):
                line = line + "  "
            for col in range (0,n-1):
                line = line + "\_"
            line = line + "\     "
            for col in range (Fork,O5):
                line = line + "\_"
            max = O3
            if (max == ns) : max = max -1
            for col in range (O5,max):
                line = line + "|_"
            for col in range (O3,ns-1):
                line = line + "/_"
            if (O3 == ns) : line = line + "|"
            else : line = line + "/"
        if (lin == 4) :
            line = line + " "
            for col in range (n,ns):
                line = line + "  "
            for col in range (0,n-1):
                line = line + "  "
            if Fork == 0 : line = line + " \   /"
            else : line = line + " \_ _/"
        if (lin > 4) :
            line = line + " "
            for col in range (n,ns):
                line = line + "  "
            for col in range (0,n-1):
                line = line + "  "
            if Fork < lin-3 : line = line + "      "
            else : line = line + " |_ _|"
        print(line)
    return

def mk_fork_lin(ns,lin,O5,Fork,O3):
    n = ns - Fork
    if O5 < Fork : O5 = Fork
    if O5 > ns-1 : O5 = ns-1
    if O3 < O5 : O3 = O5
    if O3 > ns : O3 = ns
    size = 7 + 4*ns 
    #print("Fork = %d; O5 = %d; O3 = %d" %(Fork,O5,O3))
    line = (" ")
    if (lin == 0) :
        for col in range (0,5+(2*n)):
            line = line + " "
        if (O5 > 0):
            line = line + " "
            for col in range (0,O5-1):
                line = line + "__"
            for col in range (O5-1,O5):
                line = line + "  "
        for col in range (2*O5,2*O3-1):
            line = line + " "
        if (O3 < n-1):
            line = line + " _"
            for col in range (O3+1,n-1):
                line = line + "__"
            line = line + "_ "
        l = len(line)
        for col in range (l,size):
            line = line + " "
    if (lin == 1) :
        for col in range (0,3+(2*n)):
            line = line + " "
        line = line + "5'"
        if (O5 > 0):
            for col in range (0,O5-1):
                line = line + "/ "
            for col in range (O5-1,O5):
                line = line + "/\\"
        for col in range (2*O5,2*O3-1):
            line = line + "_"
        if (O3 < ns):
            line = line + "/\\"
            for col in range (O3+1,ns):
                line = line + " \\"
        line = line + "3'"
        l = len(line)
        for col in range (l,size):
            line = line + " "
    if (lin == 2) :
        for col in range (0,5+(2*n)):
            line = line + " "
        for col in range (0,O5):
            line = line + "  "
        for col in range (O5,O3):
            line = line + "| "
        for col in range (O3,n+1):
            line = line + " "
        l = len(line)
        for col in range (l,size):
            line = line + " "
    if (lin == 3) :
        line = line + " "
        for col in range (n,ns):
            line = line + "  "
        for col in range (0,n-1):
            line = line + "\_"
        line = line + "\     "
        for col in range (Fork,O5):
            line = line + "\_"
        max = O3
        if (max == ns) : max = max -1
        for col in range (O5,max):
            line = line + "|_"
        for col in range (O3,ns-1):
            line = line + "/_"
        if (O3 == ns) : line = line + "|"
        else : line = line + "/"
        l = len(line)
        for col in range (l,size):
            line = line + " "
    if (lin == 4) :
        line = line + " "
        for col in range (n,ns):
            line = line + "  "
        for col in range (0,n-1):
            line = line + "  "
        if Fork == 0 : line = line + " \   /"
        else : line = line + " \_ _/"
        l = len(line)
        for col in range (l,size):
            line = line + " "
    if (lin > 4) :
        line = line + " "
        for col in range (n,ns):
            line = line + "  "
        for col in range (0,n-1):
            line = line + "  "
        if Fork < lin-3 : line = line + "      "
        else : line = line + " |_ _|"
        l = len(line)
        for col in range (l,size):
            line = line + " "
    return line


def disp_fork_comp(ns,O5a,Forka,O3a,O5b,Forkb,O3b):
    n = ns - Forka
    if O5a < Forka : O5a = Forka
    if O5a > ns-1 : O5a = ns-1
    if O3a < O5a : O3a = O5a
    if O3a > ns : O3a = ns
    if O5b < Forkb : O5b = Forkb
    if O5b > ns-1 : O5b = ns-1
    if O3b < O5b : O3b = O5b
    if O3b > ns : O3b = ns 
    print("Fork = %d; O5 = %d; O3 = %d // Fork = %d; O5 = %d; O3 = %d" %(Forka,O5a,O3a,Forkb,O5b,O3b))
    for lin in range (0,ns+4):
        line = mk_fork_lin(ns,lin,O5a,Forka,O3a)
        line = line + " "
        line = line + mk_fork_lin(ns,lin,O5b,Forkb,O3b)
        print(line)
    return

def disp_three_fork(ns,O5a,Forka,O3a,O5b,Forkb,O3b,O5c,Forkc,O3c):
    n = ns - Forka
    sec = third = 1;
    if O5a < Forka : O5a = Forka
    if O5a > ns-1 : O5a = ns-1
    if O3a < O5a : O3a = O5a
    if O3a > ns : O3a = ns
    if (Forkb < 0) or (O3b < 0) or (O5b < 0) : sec = 0
    if O5b < Forkb : O5b = Forkb
    if O5b > ns-1 : O5b = ns-1
    if O3b < O5b : O3b = O5b
    if O3b > ns : O3b = ns
    if (Forkc < 0) or (O3c < 0) or (O5c < 0) : third = 0
    if O5c < Forkc : O5c = Forkc
    if O5c > ns-1 : O5c = ns-1
    if O3c < O5c : O3c = O5c
    if O3c > ns : O3c = ns
    line = ("Fork = %d; O5 = %d; O3 = %d //" %(Forka,O5a,O3a))
    if sec > 0 : line = line + (" Fork = %d; O5 = %d; O3 = %d //" %(Forkb,O5b,O3b))
    if third > 0 : line = line + (" Fork = %d; O5 = %d; O3 = %d" %(Forkc,O5c,O3c))
    print(line)
    for lin in range (0,ns+4):
        line = mk_fork_lin(ns,lin,O5a,Forka,O3a)
        if sec > 0 : line = line + " " + mk_fork_lin(ns,lin,O5b,Forkb,O3b)
        if third > 0 : line = line + " " + mk_fork_lin(ns,lin,O5c,Forkc,O3c)
        print(line)
    return


def disp_oligo(n,O5,O3):
    if O5 < 0 : O5 = 0
    if O5 > n-1 : O5 = n-1
    if O3 < 0 : O3 = 0
    if O3 > n : O3 = n
    print("O5 = %d; O3 = %d" % (O5,O3))
    for lin in range (0,4):
        #line = ("%d %d " % (lin,n))
        line = (" ")
        if (lin == 0) :
            line = line + "  "
            if (O5 > 0):
                line = line + " "
                for col in range (0,O5-1):
                    line = line + "__"
                for col in range (O5-1,O5):
                    line = line + " "
            for col in range (2*O5,2*O3-1):
                line = line + " "
            if (O3 < n-1):
                line = line + " _"
                for col in range (O3+1,n-1):
                    line = line + "__"
                line = line + "_ "
        if (lin == 1) :
            line = line + "5'"
            if (O5 > 0):
                for col in range (0,O5-1):
                    line = line + "/ "
                for col in range (O5-1,O5):
                    line = line + "/\\"
            for col in range (2*O5,2*O3-1):
                line = line + "_"
            if (O3 < n):
                line = line + "/\\"
                for col in range (O3+1,n):
                    line = line + " \\"
            line = line + "3'"
        if (lin == 2) :
            line = line + "  "
            for col in range (0,O5):
                line = line + "  "
            for col in range (O5,O3):
                line = line + "| "
            for col in range (O3,n+1):
                line = line + "  "
        if (lin == 3) :
            line = line + "  "
            for col in range (0,O5):
                line = line + "\_"
            max = O3
            if (max == n) : max = max -1
            for col in range (O5,max):
                line = line + "|_"
            for col in range (O3,n-1):
                line = line + "/_"
            if (O3 == n) : line = line + "|"
            else : line = line + "/"
        print(line)
    return

def mk_oligo_lin(n,lin,O5,O3):
    if O5 < 0 : O5 = 0
    if O5 > n-1 : O5 = n-1
    if O3 < 0 : O3 = 0
    if O3 > n : O3 = n
    size = 2*(n+3)
    #print("O5 = %d; O3 = %d" % (O5,O3))
    line = (" ")
    if (lin == 0) :
        line = line + "  "
        if (O5 > 0):
            line = line + " "
            for col in range (0,O5-1):
                line = line + "__"
            for col in range (O5-1,O5):
                line = line + " "
        for col in range (2*O5,2*O3-1):
            line = line + " "
        if (O3 < n-1):
            line = line + " _"
            for col in range (O3+1,n-1):
                line = line + "__"
            line = line + "_ "
    if (lin == 1) :
        line = line + "5'"
        if (O5 > 0):
            for col in range (0,O5-1):
                line = line + "/ "
            for col in range (O5-1,O5):
                line = line + "/\\"
        for col in range (2*O5,2*O3-1):
            line = line + "_"
        if (O3 < n):
            line = line + "/\\"
            for col in range (O3+1,n):
                line = line + " \\"
            line = line + "3'"
    if (lin == 2) :
        line = line + "  "
        for col in range (0,O5):
            line = line + "  "
        for col in range (O5,O3):
            line = line + "| "
        for col in range (O3,n+1):
            line = line + "  "
    if (lin == 3) :
        line = line + "  "
        for col in range (0,O5):
            line = line + "\_"
        max = O3
        if (max == n) : max = max -1
        for col in range (O5,max):
            line = line + "|_"
        for col in range (O3,n-1):
            line = line + "/_"
        if (O3 == n) : line = line + "|"
        else : line = line + "/"
    l = len(line)
    for col in range (l,size):
        line = line + " "
    return(line)


def disp_oligo_comp(ns,O5a,O3a,O5b,O3b):
    n = ns
    if O5a < 0 : O5a = 0
    if O5a > ns-1 : O5a = ns-1
    if O3a < O5a : O3a = O5a
    if O3a > ns : O3a = ns
    if O5b < 0 : O5b = 0
    if O5b > ns-1 : O5b = ns-1
    if O3b < O5b : O3b = O5b
    if O3b > ns : O3b = ns 
    print("O5 = %d; O3 = %d // O5 = %d; O3 = %d" %(O5a,O3a,O5b,O3b))
    for lin in range (0,ns+4):
        line = mk_oligo_lin(ns,lin,O5a,O3a)
        line = line + " "
        line = line + mk_oligo_lin(ns,lin,O5b,O3b)
        print(line)
    return

def disp_three_oligo(ns,O5a,O3a,O5b,O3b,O5c,O3c):
    n = ns
    sec = third = 1;
    if O5a < 0 : O5a = 0
    if O5a > ns-1 : O5a = ns-1
    if O3a < O5a : O3a = O5a
    if O3a > ns : O3a = ns
    if (O3b < 0) or (O5b < 0) : sec = 0
    if O5b < 0 : O5b = 0
    if O5b > ns-1 : O5b = ns-1
    if O3b < O5b : O3b = O5b
    if O3b > ns : O3b = ns
    if (O3c < 0) or (O5c < 0) : third = 0
    if O5c < 0 : O5c = 0
    if O5c > ns-1 : O5c = ns-1
    if O3c < O5c : O3c = O5c
    if O3c > ns : O3c = ns
    line = ("O5 = %d; O3 = %d //" % (O5a,O3a))
    if sec > 0 : line = line + ("O5 = %d; O3 = %d " %(O5b,O3b))
    if third > 0 : line = line +  ("// O5 = %d; O3 = %d" %(O5c,O3c))
    for lin in range (0,ns+4):
        line = mk_oligo_lin(ns,lin,O5a,O3a)
        if sec > 0 : line = line + " " + mk_oligo_lin(ns,lin,O5b,O3b)
        if third > 0 : line = line + " " + mk_oligo_lin(ns,lin,O5c,O3c)
        print(line)
    return


def  fill_loop_mat(mat,ind,l,nn,rco,roo,roo2,ren,disp_flag,LAST_2bp):
#  \
#   \___/
#  ------
# 5'->i = 2  3'->j=5
# 5' end transitions possible state
    for j in range (2,nn):
        for i in range(0,j-1):
            in0 = ind[i,j];
            in1 = ind[i+1,j];
            bp_remain = j - i;
            mat[in0,in1] += rco;    # closing
            mat[in1,in1] += -rco;    # closing
            if (LAST_2bp > 0) and (bp_remain <= 2) :
                mat[in1,in0] += roo2[i];    # openning reduced energy
                mat[in0,in0] += -roo2[i];    # openning reduced energy
                if ((disp_flag > 3) and (nn < 10)):
                    print("Peeling of the oligo 5'end  on 1 base [%d,%d]->%d and [%d,%d]->%d dg lower %g bpremain %d" % (i,j,in0,i+1,j,in1,roo2[i],bp_remain));
                    disp_oligo_comp(nn-1,i,j,i+1,j)
            else :
                mat[in1,in0] += roo[i];    # openning
                mat[in0,in0] += -roo[i];    # openning
                if ((disp_flag > 3) and (nn < 10)):
                    print("Peeling of the oligo 5'end  on 1 base [%d,%d]->%d and [%d,%d]->%d dg %g bpremain %d" % (i,j,in0,i+1,j,in1,roo[i],bp_remain));
                    disp_oligo_comp(nn-1,i,j,i+1,j)


    if ((disp_flag > 3) and (nn < 8)):
        print(" \n5'end transitions possible state");
        print(mat);
        print("\n");


# 3' end transitions possible state

    for i in range(0,nn):
        for j in range (i+2,nn):
            in0 = ind[i,j];
            in1 = ind[i,j-1];
            bp_remain = j - i;
            mat[in0,in1] += rco;    # closing
            mat[in1,in1] += -rco;    # closing
            if (LAST_2bp > 0) and (bp_remain <= 2) :
                mat[in1,in0] += roo2[j-1];    # openning reduced energy
                mat[in0,in0] += -roo2[j-1];    # openning reduced energy
                if ((disp_flag > 3) and (nn < 10)):
                    print("Peeling of the oligo 3'end  on 1 base [%d,%d]->%d and [%d,%d]->%d dg lower %g bpremain %d" % (i,j,in0,i,j-1,in1,roo2[j-1],bp_remain));
                    disp_oligo_comp(nn-1,i,j,i,j-1)
            else :
                mat[in1,in0] += roo[j-1];    # openning
                mat[in0,in0] += -roo[j-1];    # openning
                if ((disp_flag > 3) and (nn < 10)):
                    print("Peeling of the oligo 3'end  on 1 base [%d,%d]->%d and [%d,%d]->%d dg %g bpremain %d" % (i,j,in0,i,j-1,in1,roo[j-1],bp_remain));
                    disp_oligo_comp(nn-1,i,j,i,j-1)

    if ((disp_flag > 3) and (nn < 8)):
        print(" \n3'end transitions possible state");
        print(mat);
        print("\n");

    for i in range(0,l):
        tmp = 0;
        for j in range(0,l):
            tmp += mat[j,i];
        if ((disp_flag > 3) and (nn < 8)):
            print ("col %d sum %g" % (i,tmp));


#escaping terms
    for i in range (0,nn-1):
        in0 = ind[i,i+1];
        if LAST_2bp > 0:
            tmp = -roo2[i];
        else :
            tmp = -roo[i];
        mat[in0,in0] += tmp;
        if ((disp_flag > 3) and (nn < 10)):
            print("5' transition escaping [%d,%d]->%d by %g" % (i,i+1,in0,tmp));
            #we need to escape also for 3' front which repeats the same process
            #    for i in range (0,nn-1):
            #        in0 = ind[i,i+1];
            #        mat[in0,in0] += -roo2[i];
            #        if ((disp_flag > 1) and (nn < 8)):
            #            print("3' transition escaping [%d,%d]->%d" % (i,i+1,in0));
    if (disp_flag > 3) and (nn < 10):
        print(" \ndiagonal terms");
        print(mat);
        print("\nEncercling transition ren = %g\n" % (ren));
        #encircling transitions
    if nn-1 > 6: 
        nenc = 3
    else :
        nenc = nn - 1
    for i in range(0,nenc): #was 3, (int)((nn-1)/2)):
        in0=ind[i,nn-i-1]
        mat[in0,in0]+=-(ren)**(nn-i)
        if (disp_flag > 2):
            print("encercling escaping [%d,%d]->%d ren %g rate %g" % (i,nn-i-1,in0,ren,(ren)**(nn-i)));
                # fork displacement
    return mat

def   fill_fork_mat(mat,ind,l,nn,rco,rch,roo,roo2,roh,disp_flag,LAST_2bp):
# 5' end transitions possible state
    for j in range (0,nn-1):
        for i in range(j,nn-2):
            for k in range (i+2,nn):
                in0 = ind[i,j,k];
                in1 = ind[i+1,j,k];
                bp_remain = k - i;
                mat[in0,in1] += rco;    # closing
                mat[in1,in1] += -rco;    # closing
                if (LAST_2bp > 0) and (bp_remain <= 2) :
                    mat[in1,in0] += roo2[i];    # openning reduced energy
                    mat[in0,in0] += -roo2[i];    # openning reduced energy
                    if ((disp_flag > 3) and (nn < 8)):
                        print("Peeling one base on 5' end of oligo [%d,%d,%d]->%d %g and [%d,%d,%d]->%d lower %g bpremain %d " % (i,j,k,in0,rco,i+1,j,k,in1,roo2[i],bp_remain));
                        disp_fork_comp(nn-1,i,j,k,i+1,j,k)
                else :
                    mat[in1,in0] += roo[i];    # openning
                    mat[in0,in0] += -roo[i];    # openning
                    if ((disp_flag > 3) and (nn < 8)):
                        print("Peeling one base on 5' end of oligo [%d,%d,%d]->%d %g and [%d,%d,%d]->%d %g bpremain %d " % (i,j,k,in0,rco,i+1,j,k,in1,roo[i],bp_remain));
                        disp_fork_comp(nn-1,i,j,k,i+1,j,k)


    if ((disp_flag > 3) and (nn < 8)):
        print(" \n5'end transitions possible state");
        print(mat);
        print("\n");


# hairpin transitions possible state
    for i in range(1,nn):
        for j in range (0,i):
            for k in range (i+1,nn):
                in0 = ind[i,j,k];
                in1 = ind[i,j+1,k];
                mat[in0,in1] += roh[j];    # openning
                mat[in1,in1] += -roh[j];    # openning
                mat[in1,in0] += rch;    # closing
                mat[in0,in0] += -rch;    # closing
                if ((disp_flag > 3) and (nn < 8)):
                    print("fork closing one base [%d,%d,%d]->%d %g and [%d,%d,%d]->%d %g" % (i,j,k,in0,rch,i,j+1,k,in1,roh[j]));
                    disp_fork_comp(nn-1,i,j,k,i,j+1,k)
    if ((disp_flag > 3) and (nn < 8)):
        print(" \nfork transitions possible state");
        print(mat);
        print("\n");

# 3' end transitions possible state
    for i in range(0,nn):
        for j in range (0,i+1):
            for k in range (i+2,nn):
                in0 = ind[i,j,k];
                in1 = ind[i,j,k-1];
                bp_remain = k - i;
                mat[in0,in1] += rco;    # closing
                mat[in1,in1] += -rco;    # closing
                if (LAST_2bp > 0) and (bp_remain <= 2) :
                    mat[in1,in0] += roo2[k-1];    # openning
                    mat[in0,in0] += -roo2[k-1];    # openning
                    if ((disp_flag > 3) and (nn < 8)):
                        print("Peeling one base on the 3'end  [%d,%d,%d]->%d and [%d,%d,%d]->%d dg lower %g bpremain %d" % (i,j,k,in0,i,j,k-1,in1,roo2[k-1],bp_remain));
                        disp_fork_comp(nn-1,i,j,k,i,j,k-1)
                else :
                    mat[in1,in0] += roo[k-1];    # openning
                    mat[in0,in0] += -roo[k-1];    # openning
                    if ((disp_flag > 3) and (nn < 8)):
                        print("Peeling one base on the 3'end [%d,%d,%d]->%d and [%d,%d,%d]->%d dg %g bpremain %d" % (i,j,k,in0,i,j,k-1,in1,roo[k-1],bp_remain));
                        disp_fork_comp(nn-1,i,j,k,i,j,k-1)

    if ((disp_flag > 3) and (nn < 8)):
        print(" \n3'end transitions possible state");
        print(mat);
        print("\n");
        for i in range(0,l):
            tmp = 0;
            for j in range(0,l):
                tmp += mat[j,i];
            print ("col %d sum %g" % (i,tmp));


#escaping terms
    for i in range (0,nn-1):
        for j in range (0,i+1):
            in0 = ind[i,j,i+1]; #  bp_remain = k - i = 1;
            if LAST_2bp > 0:
                mat[in0,in0] += -roo2[i];
            else :
                mat[in0,in0] += -roo[i];
            if ((disp_flag > 3) and (nn < 8)):
                if LAST_2bp > 0:
                    print("Escaping by opening the last base [%d,%d,%d]->%d %g" % (i,j,i+1,in0,roo2[i]));
                else :
                    print("Escaping by opening the last base [%d,%d,%d]->%d %g" % (i,j,i+1,in0,roo[i]));
                print("\n");
                disp_fork(nn-1,i,j,i+1)



#    for i in range (0,nn-1):
#        for j in range (0,i+1):
#            in0 = ind[i,j,i+1];
#            mat[in0,in0] += -roo2[i];
#            if ((disp_flag > 2) and (nn < 8)):
#                print("3' transition escaping [%d,%d,%d]->%d %g" % (i,j,i+1,in0,roo2[i]));# roo[i


    if (disp_flag > 3) and (nn < 8):
        print(" \ndiagonal terms");
        print(mat);
        print("\n");

    return mat

complementary_base = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
def is_self_complementary(seq):
    '''
    Returns true if seq is a self-complementary DNA sequence, and returns false otherwise.
    '''
    length = len(seq)
    # assumes that a sequence is self-complementary, and checks wheter that's actually the case.
    self_complementary = True
    for i in range(length):
        i_c = length - 1 - i
        if seq[i] != complementary_base[seq[i_c]]:
            self_complementary = False
            break
 
    return self_complementary




def salt_correction2(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, method=1, seq=None):
    """Calculate a term to correct Tm for salt ions.

    Depending on the Tm calculation, the term will correct Tm or entropy. To
    calculate corrected Tm values, different operations need to be applied:

     - methods 1-4: Tm(new) = Tm(old) + corr
     - method 5: deltaH(new) = deltaH(old) + corr
     - methods 6+7: Tm(new) = 1/(1/Tm(old) + corr)

    Parameters:
     - Na, K, Tris, Mg, dNTPS: Millimolar concentration of respective ion. To
       have a simple 'salt correction', just pass Na. If any of K, Tris, Mg and
       dNTPS is non-zero, a 'sodium-equivalent' concentration is calculated
       according to von Ahsen et al. (2001, Clin Chem 47: 1956-1961):
       [Na_eq] = [Na+] + [K+] + [Tris]/2 + 120*([Mg2+] - [dNTPs])^0.5
       If [dNTPs] >= [Mg2+]: [Na_eq] = [Na+] + [K+] + [Tris]/2
     - method: Which method to be applied. Methods 1-4 correct Tm, method 5
       corrects deltaS, methods 6 and 7 correct 1/Tm. The methods are:

       1. 16.6 x log[Na+]
          (Schildkraut & Lifson (1965), Biopolymers 3: 195-208)
       2. 16.6 x log([Na+]/(1.0 + 0.7*[Na+]))
          (Wetmur (1991), Crit Rev Biochem Mol Biol 126: 227-259)
       3. 12.5 x log(Na+]
          (SantaLucia et al. (1996), Biochemistry 35: 3555-3562
       4. 11.7 x log[Na+]
          (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465
       5. Correction for deltaS: 0.368 x (N-1) x ln[Na+]
          (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465)
       6. (4.29(%GC)-3.95)x1e-5 x ln[Na+] + 9.40e-6 x ln[Na+]^2
          (Owczarzy et al. (2004), Biochemistry 43: 3537-3554)
       7. Complex formula with decision tree and 7 empirical constants.
          Mg2+ is corrected for dNTPs binding (if present)
          (Owczarzy et al. (2008), Biochemistry 47: 5336-5353)

    Examples:
        >>> from Bio.SeqUtils import MeltingTemp as mt
        >>> print('%0.2f' % mt.salt_correction(Na=50, method=1))
        -21.60
        >>> print('%0.2f' % mt.salt_correction(Na=50, method=2))
        -21.85
        >>> print('%0.2f' % mt.salt_correction(Na=100, Tris=20, method=2))
        -16.45
        >>> print('%0.2f' % mt.salt_correction(Na=100, Tris=20, Mg=1.5,
        ...                                    method=2))
        -10.99

    """
    if method in (5, 6, 7, 8) and not seq:
        raise ValueError('sequence is missing (is needed to calculate ' +
                         'GC content or sequence length).')
    if seq:
        seq = str(seq)
    corr = 0
    if not method:
        return corr
    Mon = Na + K + Tris / 2.0  # Note: all these values are millimolar
    mg = Mg * 1e-3             # Lowercase ions (mg, mon, dntps) are molar
    # Na equivalent according to von Ahsen et al. (2001):
    #print("1) Mon %g mg %g" % (Mon,mg))
    if sum((K, Mg, Tris, dNTPs)) > 0 and not method == 7 and dNTPs < Mg:
        # dNTPs bind Mg2+ strongly. If [dNTPs] is larger or equal than
        # [Mg2+], free Mg2+ is considered not to be relevant.
        Mon += 120 * math.sqrt(Mg - dNTPs)
        #print("2) Mon %g mg %g" % (Mon,mg))
    mon = Mon * 1e-3
    #print("3) Mon %g mg %g" % (mon,mg))
    # Note: math.log = ln(), math.log10 = log()
    if method in range(1, 8) and not mon:
        raise ValueError('Total ion concentration of zero is not allowed in ' +
                         'this method.')
    if method == 1:
        corr = 16.6 * math.log10(mon)
    if method == 2:
        corr = 16.6 * math.log10((mon) / (1.0 + 0.7 * (mon)))
    if method == 3:
        corr = 12.5 * math.log10(mon)
    if method == 4:
        corr = 11.7 * math.log10(mon)
    if method == 5:
        corr = 0.368 * (len(seq) - 1) * math.log(mon)
    if method == 6:
        corr = (4.29 * GC(seq) / 100 - 3.95) * 1e-5 * math.log(mon) +\
            9.40e-6 * math.log(mon) ** 2
    if method == 7:
        a, b, c, d = 3.92, -0.911, 6.26, 1.42
        e, f, g = -48.2, 52.5, 8.31
        if dNTPs > 0:
            dntps = dNTPs * 1e-3
            ka = 3e4  # Dissociation constant for Mg:dNTP
            # Free Mg2+ calculation:
            mg = (-(ka * dntps - ka * mg + 1.0) +
                  math.sqrt((ka * dntps - ka * mg + 1.0) ** 2 +
                            4.0 * ka * mg)) / (2.0 * ka)
        if Mon > 0:
            R = math.sqrt(mg) / mon
            if R < 0.22:
                corr = (4.29 * SeqUtils.GC(seq) / 100 - 3.95) * \
                    1e-5 * math.log(mon) + 9.40e-6 * math.log(mon) ** 2
                return corr
            elif R < 6.0:
                a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
                d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) -
                            8.03e-3 * math.log(mon) ** 2)
                g = 8.31 * (0.486 - 0.258 * math.log(mon) +
                            5.25e-3 * math.log(mon) ** 3)
        corr = (a + b * math.log(mg) + (SeqUtils.GC(seq) / 100) *
                (c + d * math.log(mg)) + (1 / (2.0 * (len(seq) - 1))) *
                (e + f * math.log(mg) + g * math.log(mg) ** 2)) * 1e-5
    if method == 8:  # to compare to santalucia.py
        corr = 0.368 * len(seq) * math.log(mon)
    if method > 8:
        raise ValueError('Allowed values for parameter \'method\' are 1-7.')
    #print("salt corr %g" % (corr))
    return corr

# SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440, corrigé VC
DNA_DNA_NN4 = {
    "init": (0.2, -5.7), "init_A/T": (2.2, 6.9), "init_G/C": (0, 0),
    "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
    "sym": (0, -1.4),
    "AA/TT": (-7.6, -21.3), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9)}


# RNA/DNA
# Sugimoto et al. (1995), Biochemistry 34: 11211-11216
RNA_DNA_NN1 = {
    'init': (1.9, -3.9), 'init_A/T': (0, 0), 'init_G/C': (0, 0),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, 0),
    'AA/TT': (-11.5, -36.4), 'AC/TG': (-7.8, -21.6), 'AG/TC': (-7.0, -19.7),
    'AT/TA': (-8.3, -23.9), 'CA/GT': (-10.4, -28.4), 'CC/GG': (-12.8, -31.9),
    'CG/GC': (-16.3, -47.1), 'CT/GA': (-9.1, -23.5), 'GA/CT': (-8.6, -22.9),
    'GC/CG': (-8.0, -17.1), 'GG/CC': (-9.3, -23.2), 'GT/CA': (-5.9, -12.3),
    'TA/AT': (-7.8, -23.2), 'TC/AG': (-5.5, -13.5), 'TG/AC': (-9.0, -26.1),
    'TT/AA': (-7.8, -21.9)}


# LNA/DNA
# Owczarzy R, You Y, Groth CL, and Tataurov AV (2011) Stability and mismatch discrimination of locked nucleic acid-DNA duplexes. Biochemistry, 50:9352-9367.
LNA_DNA_NN5 = {
    'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),
    'AA/TT':(-9.991,-27.175),'AC/TG':(-11.389,-28.963),'AG/TC':(-12.793,-31.607),'AT/TA':(-14.703,-40.750),
    'CA/GT':(-14.177,-35.498),'CC/GG':(-15.399,-36.375),'CG/GC':(-14.558,-35.239),'CT/GA':(-15.737,-41.218),
    'GA/CT':(-13.959,-35.097),'GC/CG':(-16.109,-40.738),'GG/CC':(-13.022,-29.673),'GT/CA':(-17.361,-45.858),
    'TA/AT':(-10.318,-26.108),'TC/AG':(-9.166,-21.535),'TG/AC':(-10.046,-22.591),'TT/AA':(-10.419,-27.683)}






LNA_DNA_NN1 = {
    'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),

    #'aa/tt':(-7.9,-22.25),'ac/tg':(-8.4,-22.45),'ag/tc':(-7.8,-21.025),'at/ta':(-7.2,-20.375),
    #'ca/gt':(-8.5,-22.725),'cc/gg':(-8.0,-19.85),'cg/gc':(-10.6,-27.2),'ct/ga':(-7.8,-21.025),
    #'ga/ct':(-8.2,-22.25),'gc/cg':(-9.8,-24.375),'gg/cc':(-8.0,-19.85),'gt/ca':(-8.4,-22.45),
    #'ta/at':(-7.2,-21.35),'tc/ag':(-8.2,-22.25),'tg/ac':(-8.5,-22.725),'tt/aa':(-7.9,-22.25),

    'aa/tt':(-7.9,-22.2),'ac/tg':(-8.4,-22.4),'ag/tc':(-7.8,-21.0),'at/ta':(-7.2,-20.4),
    'ca/gt':(-8.5,-22.7),'cc/gg':(-8.0,-19.9),'cg/gc':(-10.6,-27.2),'ct/ga':(-7.8,-21.0),
    'ga/ct':(-8.2,-22.2),'gc/cg':(-9.8,-24.4),'gg/cc':(-8.0,-19.9),'gt/ca':(-8.4,-22.4),
    'ta/at':(-7.2,-21.3),'tc/ag':(-8.2,-22.2),'tg/ac':(-8.5,-22.7),'tt/aa':(-7.9,-22.2),

    'AA/tt':(-9.991,-27.175),'AC/tg':(-11.389,-28.963),'AG/tc':(-12.793,-31.607),'AT/ta':(-14.703,-40.750),
    'CA/gt':(-14.177,-35.498),'CC/gg':(-15.399,-36.375),'CG/gc':(-14.558,-35.239),'CT/ga':(-15.737,-41.218),
    'GA/ct':(-13.959,-35.097),'GC/cg':(-16.109,-40.738),'GG/cc':(-13.022,-29.673),'GT/ca':(-17.361,-45.858),
    'TA/at':(-10.318,-26.108),'TC/ag':(-9.166,-21.535),'TG/ac':(-10.046,-22.591),'TT/aa':(-10.419,-27.683),


    'Aa/tt':(-7.193,-19.723),'Ac/tg':(-7.269,-18.336),'Ag/tc':(-7.536,-18.387),'At/ta':(-4.918,-12.943),
    'Ca/gt':(-7.451,-18.380),'Cc/gg':(-5.904,-11.904),'Cg/gc':(-9.815,-23.491),'Ct/ga':(-7.092,-16.825),
    'Ga/ct':(-5.038,11.656),'Gc/cg':(-10.160,-24.651),'Gg/cc':(-10.844,-26.580),'Gt/ca':(-8.612,-22.327),
    'Ta/at':(-7.246,-19.738),'Tc/ag':(-6.307,-15.515),'Tg/ac':(-10.040,-25.744),'Tt/aa':(-6.372,-16.902),

    'aA/tt':(-6.908,-18.135),'aC/tg':(-5.510,-11.824),'aG/tc':(-9.000,-22.826),'aT/ta':(-5.384,-13.537),
    'cA/gt':(-7.142,-18.333),'cC/gg':(-5.937,-12.335),'cG/gc':(-10.876,-27.918),'cT/ga':(-9.471,-25.070),
    'gA/ct':(-7.756,-19.302),'gC/cg':(-10.725,-25.511),'gG/cc':(-8.943,-20.833),'gT/ca':(-9.035,-22.742),
    'tA/at':(-5.609,-16.019),'tC/ag':(-7.591,-19.031),'tG/ac':(-6.335,-15.537),'tT/aa':(-5.574,-14.149)}

LNA_DNA_NN2 = {
    'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),

    #'aa/tt':(-7.9,-22.25),'ac/tg':(-8.4,-22.45),'ag/tc':(-7.8,-21.025),'at/ta':(-7.2,-20.375),
    #'ca/gt':(-8.5,-22.725),'cc/gg':(-8.0,-19.85),'cg/gc':(-10.6,-27.2),'ct/ga':(-7.8,-21.025),
    #'ga/ct':(-8.2,-22.25),'gc/cg':(-9.8,-24.375),'gg/cc':(-8.0,-19.85),'gt/ca':(-8.4,-22.45),
    #'ta/at':(-7.2,-21.35),'tc/ag':(-8.2,-22.25),'tg/ac':(-8.5,-22.725),'tt/aa':(-7.9,-22.25),

    'aa/tt':(-7.9,-22.2),'ac/tg':(-8.4,-22.4),'ag/tc':(-7.8,-21.0),'at/ta':(-7.2,-20.4),
    'ca/gt':(-8.5,-22.7),'cc/gg':(-8.0,-19.9),'cg/gc':(-10.6,-27.2),'ct/ga':(-7.8,-21.0),
    'ga/ct':(-8.2,-22.2),'gc/cg':(-9.8,-24.4),'gg/cc':(-8.0,-19.9),'gt/ca':(-8.4,-22.4),
    'ta/at':(-7.2,-21.3),'tc/ag':(-8.2,-22.2),'tg/ac':(-8.5,-22.7),'tt/aa':(-7.9,-22.2),

    'AA/tt':(-9.991,-27.175),'AC/tg':(-11.389,-28.963),'AG/tc':(-12.793,-31.607),'AT/ta':(-14.703,-40.750),
    'CA/gt':(-14.177,-35.498),'CC/gg':(-15.399,-36.375),'CG/gc':(-14.558,-35.239),'CT/ga':(-15.737,-41.218),
    'GA/ct':(-13.959,-35.097),'GC/cg':(-16.109,-40.738),'GG/cc':(-13.022,-29.673),'GT/ca':(-17.361,-45.858),
    'TA/at':(-10.318,-26.108),'TC/ag':(-9.166,-21.535),'TG/ac':(-10.046,-22.591),'TT/aa':(-10.419,-27.683),


    'Aa/tt':(-7.193,-19.723),'Ac/tg':(-7.269,-18.336),'Ag/tc':(-7.536,-18.387),'At/ta':(-4.918,-12.943),
    'Ca/gt':(-7.451,-18.380),'Cc/gg':(-5.904,-11.904),'Cg/gc':(-9.815,-23.491),'Ct/ga':(-7.092,-16.825),
    'Ga/ct':(-5.038,11.656),'Gc/cg':(-10.160,-24.651),'Gg/cc':(-10.844,-26.580),'Gt/ca':(-8.612,-22.327),
    'Ta/at':(-7.246,-19.738),'Tc/ag':(-6.307,-15.515),'Tg/ac':(-10.040,-25.744),'Tt/aa':(-6.372,-16.902),

    'aA/tt':(-6.908,-18.135),'aC/tg':(-5.510,-11.824),'aG/tc':(-9.000,-22.826),'aT/ta':(-5.384,-13.537),
    'cA/gt':(-7.142,-18.333),'cC/gg':(-5.937,-12.335),'cG/gc':(-10.876,-27.918),'cT/ga':(-9.471,-25.070),
    'gA/ct':(-7.756,-19.302),'gC/cg':(-10.725,-25.511),'gG/cc':(-8.943,-20.833),'gT/ca':(-9.035,-22.742),
    'tA/at':(-5.609,-16.019),'tC/ag':(-7.591,-19.031),'tG/ac':(-6.335,-15.537),'tT/aa':(-5.574,-14.149)}


# first adjustment 
LNA_DNA_NN2O = {
    'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),

    'aa/tt':(-7.9,-22.2),'ac/tg':(-8.4,-22.4),'ag/tc':(-7.8,-21.0),'at/ta':(-7.2,-20.4),
    'ca/gt':(-8.5,-22.7),'cc/gg':(-8.0,-19.9),'cg/gc':(-10.6,-27.2),'ct/ga':(-7.8,-21.0),
    'ga/ct':(-8.2,-22.2),'gc/cg':(-9.8,-24.4),'gg/cc':(-8.0,-19.9),'gt/ca':(-8.4,-22.4),
    'ta/at':(-7.2,-21.3),'tc/ag':(-8.2,-22.2),'tg/ac':(-8.5,-22.7),'tt/aa':(-7.9,-22.2),

    'AA/tt':(-9.574,-27.679), 'AC/tg':(-11.306,-30.155), 'AG/tc':(-12.275,-31.875), 'AT/ta':(-14.707,-43.539),
    'CA/gt':(-13.941,-39.849), 'CC/gg':(-15.054,-41.891), 'CG/gc':(-14.033,-35.407), 'CT/ga':(-15.554,-43.284),
    'GA/ct':(-13.749,-38.202), 'GC/cg':(-16.191,-43.778), 'GG/cc':(-12.883,-31.383), 'GT/ca':(-15.714,-44.240),
    'TA/at':(-10.281,-27.838), 'TC/ag':(-9.074,-22.187),  'TG/ac':(-9.714,-22.876),  'TT/aa':(-10.296,-28.171), 

    'Aa/tt':(-6.830,-18.670), 'Ac/tg':(-7.076,-17.743), 'Ag/tc':(-7.698,-19.004), 'At/ta':(-6.604,-17.792),
    'Ca/gt':(-7.327,-17.718), 'Cc/gg':(-6.292,-12.582), 'Cg/gc':(-10.368,-24.587), 'Ct/ga':(-8.236,-19.535),
    'Ga/ct':(-5.372,-12.581), 'Gc/cg':(-9.272,-22.585), 'Gg/cc':(-9.633,-23.522), 'Gt/ca':(-8.922,-22.854),
    'Ta/at':(-6.349,-17.232), 'Tc/ag':(-6.775,-16.606), 'Tg/ac':(-8.721,-22.252), 'Tt/aa':(-7.609,-20.234),

    'aA/tt':(-7.377,-19.474), 'aC/tg':(-6.501,-14.295), 'aG/tc':(-7.606,-19.295), 'aT/ta':(-5.519,-14.000),
    'cA/gt':(-7.949,-20.088), 'cC/gg':(-6.362,-13.462), 'cG/gc':(-10.091,-25.966), 'cT/ga':(-9.066,-24.107),
    'gA/ct':(-7.545,-18.411), 'gC/cg':(-9.207,-21.816), 'gG/cc':(-8.963,-21.053), 'gT/ca':(-8.808,-22.432),
    'tA/at':(-8.327,-23.885), 'tC/ag':(-9.565,-24.129), 'tG/ac':(-7.220,-17.648), 'tT/aa':(-7.093,-18.173)}

LNA_DNA_NNO = {
    'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),

    'aa/tt':(-7.9,-22.2),'ac/tg':(-8.4,-22.4),'ag/tc':(-7.8,-21.0),'at/ta':(-7.2,-20.4),
    'ca/gt':(-8.5,-22.7),'cc/gg':(-8.0,-19.9),'cg/gc':(-10.6,-27.2),'ct/ga':(-7.8,-21.0),
    'ga/ct':(-8.2,-22.2),'gc/cg':(-9.8,-24.4),'gg/cc':(-8.0,-19.9),'gt/ca':(-8.4,-22.4),
    'ta/at':(-7.2,-21.3),'tc/ag':(-8.2,-22.2),'tg/ac':(-8.5,-22.7),'tt/aa':(-7.9,-22.2),

    'AA/tt':(-9.991,-27.175),'AC/tg':(-11.389,-28.963),'AG/tc':(-12.793,-31.607),'AT/ta':(-14.703,-40.750),
    'CA/gt':(-14.177,-35.498),'CC/gg':(-15.399,-36.375),'CG/gc':(-14.558,-35.239),'CT/ga':(-15.737,-41.218),
    'GA/ct':(-13.959,-35.097),'GC/cg':(-16.109,-40.738),'GG/cc':(-13.022,-29.673),'GT/ca':(-17.361,-45.858),
    'TA/at':(-10.318,-26.108),'TC/ag':(-9.166,-21.535),'TG/ac':(-10.046,-22.591),'TT/aa':(-10.419,-27.683),

    'Aa/tt':(-6.726,-18.443), 'Ac/tg':(-7.052,-17.790), 'Ag/tc':(-7.738,-18.880), 'At/ta':(-6.729,-17.710),
    'Ca/gt':(-7.312,-18.038), 'Cc/gg':(-6.269,-12.641), 'Cg/gc':(-10.409,-24.913), 'Ct/ga':(-7.901,-18.745),
    'Ga/ct':(-5.471,-12.659), 'Gc/cg':(-9.272,-22.497), 'Gg/cc':(-9.580,-23.482), 'Gt/ca':(-8.898,-23.069),
    'Ta/at':(-6.338,-17.265), 'Tc/ag':(-6.795,-16.714), 'Tg/ac':(-8.671,-22.235), 'Tt/aa':(-7.610,-20.185),

    'aA/tt':(-7.379,-19.371), 'aC/tg':(-6.589,-14.140), 'aG/tc':(-7.598,-19.269), 'aT/ta':(-5.543,-13.937),
    'cA/gt':(-7.913,-20.313), 'cC/gg':(-6.533,-13.573), 'cG/gc':(-10.151,-26.058), 'cT/ga':(-9.026,-23.892),
    'gA/ct':(-7.527,-18.731), 'gC/cg':(-9.185,-21.849), 'gG/cc':(-8.989,-20.941), 'gT/ca':(-8.846,-22.267),
    'tA/at':(-8.338,-23.812), 'tC/ag':(-9.597,-24.060), 'tG/ac':(-7.189,-17.633), 'tT/aa':(-7.123,-18.081 )}

'''

'Aa/tt':(-6.726,-18.443),
'Ac/':(-7.052,-17.790),
'Ag/':(-7.738,-18.880),
'At/':(-6.729,-17.710),
'Ca/':(-7.312,-18.038),
'Cc/':(-6.269,-12.641),
'Cg/':(-10.409,-24.913),
'Ct/':(-7.901,-18.745),
'Ga/':(-5.471,-12.659),
'Gc/':(-9.272,-22.497),
'Gg/':(-9.580,-23.482),
'Gt/':(-8.898,-23.069),
'Ta/':(-6.338,-17.265),
'Tc/':(-6.795,-16.714),
'Tg/':(-8.671,-22.235),
'Tt/':(-7.610,-20.185),
'aA/':(-7.379,-19.371),
'aC/':(-6.589,-14.140),
'aG/':(-7.598,-19.269),
'aT/':(-5.543,-13.937),
'cA/':(-7.913,-20.313),
'cC/':(-6.533,-13.573),
'cG/':(-10.151,-26.058),
'cT/':(-9.026,-23.892),
'gA/':(-7.527,-18.731),
'gC/':(-9.185,-21.849),
'gG/':(-8.989,-20.941),
'gT/':(-8.846,-22.267),
'tA/':(-8.338,-23.812),
'tC/':(-9.597,-24.060),
'tG/':(-7.189,-17.633),
'tT/':(-7.123,-18.081 ),

'''
# LNA/DNA
# Owczarzy R, You Y, Groth CL, and Tataurov AV (2011) Stability and mismatch discrimination of locked nucleic acid-DNA duplexes. Biochemistry, 50:9352-9367.
LNA_DNA_IMM1 = {
    'AA/at': (-3.826, -13.109),
    'AC/ag': (-2.367, -7.322),
    'AG/ac': (-4.849, -13.007),
    'AT/aa': (-5.049, -17.514),
    'AA/ta': (-4.229, -15.160),
    'CA/ga': (-5.878, -17.663),
    'GA/ca': (-8.558, -23.976),
    'TA/aa': (2.074, 3.446),
    'CA/ct': (2.218, 4.750),
    'CC/cg': (1.127, 1.826),
    'CG/cc': (-10.903, -32.025),
    'CT/ca': (-2.053, -10.517),
    'AC/tc': (1.065, -1.403),
    'CC/gc': (-9.522, -27.024),
    'GC/cc': (-4.767, -14.897),
    'TC/ac': (4.114, 9.258),
    'GA/gt': (-2.920, -9.387),
    'GC/gg': (-8.139, -21.784),
    'GG/gc': (-5.149, -12.508),
    'GT/ga': (-8.991, -27.311),
    'AG/tg': (-4.980, -15.426),
    'CG/gg': (-4.441, -12.158),
    'GG/cg': (-13.505, -36.021),
    'TG/ag': (-2.775, -9.286),
    'TA/tt': (-3.744, -12.149),
    'TC/tg': (-4.387, -13.520),
    'TG/tc': (-6.346, -16.629),
    'TT/ta': (-7.697, -25.049),
    'AT/tt': (-4.207, -14.307),
    'CT/gt': (-8.176, -22.962),
    'GT/ct': (-7.241, -20.622),
    'TT/at': (-2.051, -7.055),
    'AA/ct': (-1.362, -5.551),
    'AC/cg': (-1.759, -6.511),
    'AG/cc': (-6.549, -18.073),
    'AT/ca': (-3.563, -14.105),
    'AA/tc': (-2.078, -10.088),
    'CA/gc': (-5.868, -16.952),
    'GA/cc': (-8.477, -24.565),
    'TA/ac': (2.690, 4.965),
    'CA/at': (-9.844, -29.673),
    'CC/ag': (-3.761, -11.204),
    'CG/ac': (-9.845, -27.316),
    'CT/aa': (-3.389, -12.517),
    'AC/ta': (0.753, -0.503),
    'CC/ga': (-12.714, -35.555),
    'GC/ca': (-12.658, -35.729),
    'TC/aa': (-1.719, -7.023),
    'AA/gt': (2.193, 4.374),
    'AC/gg': (-8.453, -22.672),
    'AG/gc': (-1.164, -2.532),
    'AT/ga': (-7.418, -24.066),
    'AA/tg': (-1.963, -9.013),
    'CA/gg': (-8.712, -23.779),
    'GA/cg': (-7.875, -21.661),
    'TA/ag': (3.207, 7.156),
    'GA/at': (-2.914, -9.402),
    'GC/ag': (-9.131, -25.347),
    'GG/ac': (-2.154, -3.871),
    'GT/aa': (-8.515, -26.313),
    'AG/ta': (-6.691, -21.148),
    'CG/ga': (-3.960, -10.588),
    'GG/ca': (-12.898, -34.656),
    'TG/aa': (0.334, -0.440),
    'CA/tt': (0.382, -0.579),
    'CC/tg': (-2.716, -8.000),
    'CG/tc': (-10.363, -29.315),
    'CT/ta': (-5.783, -20.173),
    'AC/tt': (-0.692, -5.278),
    'CC/gt': (-10.299, -28.503),
    'GC/ct': (-9.062, -26.356),
    'TC/at': (2.073, 3.968),
    'TA/ct': (-5.485, -17.347),
    'TC/cg': (1.451, 1.556),
    'TG/cc': (-7.213, -20.128),
    'TT/ca': (-2.397, -11.371),
    'AT/tc': (-0.633, -5.801),
    'CT/gc': (-6.868, -21.000),
    'GT/cc': (-5.853, -16.643),
    'TT/ac': (0.211, -1.446),
    'GA/tt': (-5.551, -15.398),
    'GC/tg': (-14.943, -40.148),
    'GG/tc': (-8.110, -18.349),
    'GT/ta': (-14.213, -40.041),
    'AG/tt': (-7.130, -20.786),
    'CG/gt': (-14.862, -39.430),
    'GG/ct': (-14.622, -37.510),
    'TG/at': (-6.703, -18.111),
    'TA/gt': (-4.612, -14.039),
    'TC/gg': (-9.798, -26.406),
    'TG/gc': (-4.519, -11.065),
    'TT/ga': (-4.523, -15.693),
    'AT/tg': (-2.364, -8.834),
    'CT/gg': (-11.396, -30.732),
    'GT/cg': (-6.233, -15.933),
    'TT/ag': (-2.960, -9.305)}

# copy of DNA terminal mismatches partially corrected...
LNA_DNA_TMM1 = {
    'aa/TA': (-3.1, -7.8), 'ta/AA': (-2.5, -6.3), 'ca/GA': (-4.3, -10.7),
    'ga/CA': (-4.66, -11.43),
    'ac/TC': (-0.1, 0.5), 'tc/AC': (2.1, 4.76), 'cc/GC': (-2.1, -5.1),
    'gc/CC': (-9.62, -24.2),
    'ag/TG': (-1.1, -2.1), 'tg/AG': (-1.1, -2.7), 'cg/GG': (-3.8, -9.5),
    'gg/CG': (-8.24, -21.54),
    'at/TT': (-2.4, -6.5), 'tt/AT': (-3.2, -8.9), 'ct/GT': (-6.1, -16.9),
    'gt/CT': (-7.4, -21.2),
    'aa/TC': (-1.6, -4.0), 'ac/TA': (-1.8, -3.8), 'ca/GC': (-2.6, -5.9),
    'cc/GA': (-2.7, -6.0), 'ga/CC': (-12.8, -32.75), 'gc/CA': (-4.062, -10.62),
    'ta/AC': (-10, -30.7), 'tc/AA': (-2.7, -7.0),
    'ac/TT': (-0.9, -1.7), 'at/TC': (-2.3, -6.3), 'cc/GT': (-3.2, -8.0),
    'ct/GC': (-3.9, -10.6), 'gc/CT': (-4.9, -13.5), 'gt/CC': (-10.4, -25.7),
    'tc/AT': (-2.5, -6.3), 'tt/AC': (0.28, -1.02),
    'aa/TG': (-1.9, -4.4), 'ag/TA': (-2.5, -5.9), 'ca/GG': (-3.9, -9.6),
    'cg/GA': (-6.0, -15.5), 'ga/CG': (-9.24, -25.57), 'gg/CA': (-10.75, -26.8),
    'ta/AG': (-2.0, -4.7), 'tg/AA': (-2.4, -5.8),
    'ag/TT': (-3.2, -8.7), 'at/TG': (-3.5, -9.4), 'cg/GT': (-3.8, -9.0),
    'ct/GG': (-6.6, -18.7), 'gg/CT': (-5.7, -15.9), 'gt/CG': (-15.05, -40.65),
    'tg/AT': (-3.9, -10.5), 'tt/AG': (-3.6, -9.8)}



def Tm_NN2(seq, check=True, strict=True, c_seq=None, shift=0, nn_table=mt.DNA_NN3,
          tmm_table=mt.DNA_TMM1, imm_table=mt.DNA_IMM1, de_table=mt.DNA_DE1,
          dnac1=25, dnac2=25, selfcomp=False, Na=50, K=0, Tris=0, Mg=0,
           dNTPs=0, saltcorr=5, T_dG=25, dsCharge=0.2, F=8.5, disp_flag=2,encerclingCoeff=1,
           loop=1, encercling=1, rate=(4/1.4), oligotype = 'DNA', encerclingBpSalt = 3.0,
           dss=0.542,bpss=2.14,Sss=216,NNEnergyMod='',LAST_2bp=1):
    def GC(sequence):
        return 100 * gc_fraction(sequence, ambiguous="remove")
    """Return the Tm using nearest neighbor thermodynamics.

    Arguments:
     - seq: The primer/probe sequence as string or Biopython sequence object.
       For RNA/DNA hybridizations seq must be the RNA sequence.
     - c_seq: Complementary sequence. The sequence of the template/target in
       3'->5' direction. c_seq is necessary for mismatch correction and
       dangling-ends correction. Both corrections will automatically be
       applied if mismatches or dangling ends are present. Default=None.
     - shift: Shift of the primer/probe sequence on the template/target
       sequence, e.g.::

                           shift=0       shift=1        shift= -1
        Primer (seq):      5' ATGC...    5'  ATGC...    5' ATGC...
        Template (c_seq):  3' TACG...    3' CTACG...    3'  ACG...

       The shift parameter is necessary to align seq and c_seq if they have
       different lengths or if they should have dangling ends. Default=0
     - table: Thermodynamic NN values, eight tables are implemented:
       For DNA/DNA hybridizations:

        - DNA_NN1: values from Breslauer et al. (1986)
        - DNA_NN2: values from Sugimoto et al. (1996)
        - DNA_NN3: values from Allawi & SantaLucia (1997) (default)
        - DNA_NN4: values from SantaLucia & Hicks (2004)

       For RNA/RNA hybridizations:

        - RNA_NN1: values from Freier et al. (1986)
        - RNA_NN2: values from Xia et al. (1998)
        - RNA_NN3: valuse from Chen et al. (2012)

       For RNA/DNA hybridizations:

        - R_DNA_NN1: values from Sugimoto et al. (1995)

       Use the module's maketable method to make a new table or to update one
       one of the implemented tables.
     - tmm_table: Thermodynamic values for terminal mismatches.
       Default: DNA_TMM1 (SantaLucia & Peyret, 2001)
     - imm_table: Thermodynamic values for internal mismatches, may include
       insosine mismatches. Default: DNA_IMM1 (Allawi & SantaLucia, 1997-1998;
       Peyret et al., 1999; Watkins & SantaLucia, 2005)
     - de_table: Thermodynamic values for dangling ends:

        - DNA_DE1: for DNA. Values from Bommarito et al. (2000). Default
        - RNA_DE1: for RNA. Values from Turner & Mathews (2010)

     - dnac1: Concentration of the higher concentrated strand [nM]. Typically
       this will be the primer (for PCR) or the probe. Default=25.
     - dnac2: Concentration of the lower concentrated strand [nM]. In PCR this
       is the template strand which concentration is typically very low and may
       be ignored (dnac2=0). In oligo/oligo hybridization experiments, dnac1
       equals dnac1. Default=25.
       MELTING and Primer3Plus use k = [Oligo(Total)]/4 by default. To mimic
       this behaviour, you have to divide [Oligo(Total)] by 2 and assign this
       concentration to dnac1 and dnac2. E.g., Total oligo concentration of
       50 nM in Primer3Plus means dnac1=25, dnac2=25.
     - selfcomp: Is the sequence self-complementary? Default=False. If 'True'
       the primer is thought binding to itself, thus dnac2 is not considered.
     - Na, K, Tris, Mg, dNTPs: See method 'Tm_GC' for details. Defaults: Na=50,
       K=0, Tris=0, Mg=0, dNTPs=0.
     - saltcorr: See method 'Tm_GC'. Default=5. 0 means no salt correction.

    """
    if disp_flag > 2:
        print("Seq: %s C_seq %s oligotype -%s-" % (seq, c_seq,oligotype))
    c_seq_defined = 1
    seq = str(seq)
    if not c_seq:
        c_seq_defined = 0
        # c_seq must be provided by user if dangling ends or mismatches should
        # be taken into account. Otherwise take perfect complement.
        c_seq = Seq.Seq(seq).complement()
    c_seq = str(c_seq)
    #hairpin complementary
    h_seq = Seq.Seq(c_seq).complement()  # the HP is always without mismatch
    h_seq = str(h_seq)


    if oligotype == 'LNA' or oligotype == 'LNA2':
        h_seq = h_seq.lower();
        c_seq = c_seq.lower();
    if oligotype == 'DNA' or oligotype == 'DNA2':
        h_seq = h_seq.upper();
        c_seq = c_seq.upper();
        seq = seq.upper();
    if disp_flag > 2:
        print("Seq: %s C_seq %s " % (seq, c_seq))
    if  check and oligotype != 'LNA' and oligotype != 'LNA2':
        seq = mt._check(seq, 'Tm_NN')
        c_seq = mt._check(c_seq, 'Tm_NN')
    if disp_flag > 2:
        print("Seq: %s C_seq %s " % (seq, c_seq))
    tmp_seq = seq
    if oligotype == 'LNA' or oligotype == 'LNA2':
        tmp_cseq = c_seq.lower()
    else :
        tmp_cseq = c_seq
    tmp_hp_seq = h_seq
    if disp_flag > 2:
        print("Seq: %s C_seq %s HP_seq %s" % (tmp_seq, tmp_cseq, tmp_hp_seq))
    nn=len(seq);
    n = nn-2
    dg = numpy.zeros(nn-1);   # energy between the oligo and template (can be LNA/DNA)
    dgh = numpy.zeros(nn-1);  # energy between the 2 HP strands (it is DNA/DNA)
    T0 = 273.15
    TG = T_dG + T0
    R = 1.9872; # SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440.
    R = 1.987 # from MeltingTemp.py
    R = 1.9858775 # from SantaLucia.py
    cor = 1000.0/(R*(TG)); # cal/K-mol
    mod_erg = ""
    energy_mod = []
    if len(NNEnergyMod) > 0:
        tmp = NNEnergyMod.split(',')
        if len(tmp) != (len(seq)-1) :
            print("The energy modifier has a different size %d than the oligonucleotide %d"
                  % (len(tmp),(len(seq)-1)))
        mod_erg = "NN modified by:"
        for i in range(len(tmp)):
            energy_mod = numpy.append(energy_mod,cor*float(tmp[i]))
            mod_erg = ("%s %s," %(mod_erg,tmp[i]))
            if disp_flag > 1:
                print("for NN %d energy modified by %g kBT" % (i,energy_mod[i]))
        mod_erg = ("%s (kBT)\n" %(mod_erg))

    if ((disp_flag > 1)):
        print("Oligo sequence         %s" % (tmp_seq))
        print("Complementary sequence %s" % (tmp_cseq))
        print("Hairpin sequnece       %s" % (tmp_hp_seq))

    off_bp = 0
    delta_h = 0
    delta_s = 0
    d_h = 0  # Names for indexes
    d_s = 1  # 0 and 1
    # Dangling ends?
    if shift or len(seq) != len(c_seq):
        # Align both sequences using the shift parameter
        if shift > 0:
            tmp_seq = '.' * shift + seq
        if shift < 0:
            tmp_cseq = '.' * abs(shift) + c_seq
        if len(tmp_cseq) > len(tmp_seq):
            tmp_seq += (len(tmp_cseq) - len(tmp_seq)) * '.'
        if len(tmp_cseq) < len(tmp_seq):
            tmp_cseq += (len(tmp_seq) - len(tmp_cseq)) * '.'
        # Remove 'over-dangling' ends
        while tmp_seq.startswith('..') or tmp_cseq.startswith('..'):
            tmp_seq = tmp_seq[1:]
            tmp_cseq = tmp_cseq[1:]
        while tmp_seq.endswith('..') or tmp_cseq.endswith('..'):
            tmp_seq = tmp_seq[:-1]
            tmp_cseq = tmp_cseq[:-1]
        # Now for the dangling ends
        if tmp_seq.startswith('.') or tmp_cseq.startswith('.'):
            left_de = tmp_seq[:2] + '/' + tmp_cseq[:2]
            try:
                delta_h += de_table[left_de][d_h]
                delta_s += de_table[left_de][d_s]
            except KeyError:
                _key_error(left_de, strict)
            tmp_seq = tmp_seq[1:]
            tmp_cseq = tmp_cseq[1:]
        if tmp_seq.endswith('.') or tmp_cseq.endswith('.'):
            right_de = tmp_cseq[-2:][::-1] + '/' + tmp_seq[-2:][::-1]
            try:
                delta_h += de_table[right_de][d_h]
                delta_s += de_table[right_de][d_s]
            except KeyError:
                _key_error(right_de, strict)
            tmp_seq = tmp_seq[:-1]
            tmp_cseq = tmp_cseq[:-1]
    if delta_h > 0 or delta_s > 0 :
        if ((disp_flag > 2)):
            print ("Dangling ends dH = %g; dS = %G;" %(delta_h,delta_s))

    # Now for terminal mismatches
    #left_tmm = tmp_cseq[:2][::-1] + '/' + tmp_seq[:2][::-1]
    left_tmm = tmp_cseq[:2][::-1] + '/' + tmp_seq[:2][::-1]
    #left_tmm = left_tmm.upper()
    modif = 0
    if left_tmm in tmm_table:
        tdh = tmm_table[left_tmm][d_h]
        tds = tmm_table[left_tmm][d_s]
        delta_h += tdh
        delta_s += tds 
        tmp_seq = tmp_seq[1:]
        tmp_cseq = tmp_cseq[1:]
        tmp_hp_seq = tmp_hp_seq[1:]
        off_bp = 1
        modif = 1
        dg[0] += -cor*(tdh - (tds * TG/1000))
        if ((disp_flag > 2)):
            print ("Terminal mismatch %s  dH += %5.3f dS += %5.3f dG += %5.3f" %
                   (left_tmm, tdh, tds, tdh - (tds * TG/1000)))
    elif disp_flag > 2:
        print ("No Terminal mismatch at %s" % left_tmm)
    if oligotype == 'LNA' or oligotype == 'LNA2' :
        right_tmm = tmp_cseq[-2:] + '/' + tmp_seq[-2:]
    else :
        right_tmm = tmp_seq[-2:] + '/' + tmp_cseq[-2:]
    #right_tmm = right_tmm.upper()
    if right_tmm in tmm_table:
        tdh = tmm_table[right_tmm][d_h]
        tds = tmm_table[right_tmm][d_s]
        delta_h += tdh
        delta_s += tds
        tmp_seq = tmp_seq[:-1]
        tmp_cseq = tmp_cseq[:-1]
        tmp_hp_seq = tmp_hp_seq[:-1]
        modif = 1
        dg[nn-2] += -cor*(tdh - (tds * TG/1000))
        if ((disp_flag > 2)):
            print ("Terminal mismatch %s  dH += %5.3f dS += %5.3f dG += %5.3f" %
                   (right_tmm, tdh, tds, tdh - (tds * TG/1000)))
    elif disp_flag > 2:
        print ("No Terminal mismatch at %s" % right_tmm)
    #print("after terminal mismatch   dg[%d] = %g kBT" % (nn-2,dg[nn-2]))

    # Now everything 'unusual' at the ends is handled and removed and we can
    # look at the initiation.
    # One or several of the following initiation types may apply:

    # Type: General initiation value
    tdh = nn_table['init'][d_h]
    tds = nn_table['init'][d_s]
    delta_h += tdh
    delta_s += tds
    if tdh != 0 or tds != 0:
        if ((disp_flag > 1)):
            print ("General initiation dH = %g, added %g kcal/mol; dS = %G, added %g cal/mol;" %(delta_h,tdh,delta_s,tds))

    # Type: Duplex with no (allA/T) or at least one (oneG/C) GC pair
    sequ = seq.upper()
    if GC(sequ) == 0.0:
        tdh = nn_table['init_allA/T'][d_h]
        tds = nn_table['init_allA/T'][d_s]
        delta_h += tdh
        delta_s += tds
        if tdh != 0 or tds != 0:
            if ((disp_flag > 2)):
                print ("Duplex with (allA/T) or not a single (oneG/C) GC pair dH = %g, added %g kcal/mol; dS = %G, added %g cal/mol;"
                       %(delta_h,tdh,delta_s,tds))
    else:
        tdh = nn_table['init_oneG/C'][d_h]
        tds = nn_table['init_oneG/C'][d_s]
        delta_h += tdh
        delta_s += tds
        if tdh != 0 or tds != 0:
            if ((disp_flag > 2)):
                print ("Duplex with no (allA/T) or at least one (oneG/C) GC pair dH = %g, added %g kcal/mol; dS = %G, added %g cal/mol;"
                       %(delta_h,tdh,delta_s,tds))

    # Type: Penalty if 5' end is T
    if sequ.startswith('T'):
        tdh = nn_table['init_5T/A'][d_h]
        tds = nn_table['init_5T/A'][d_s]
        delta_h += tdh
        delta_s += tds
        if tdh != 0 or tds != 0:
            if ((disp_flag > 2)):
                print ("Penalty if 5' end is T dH = %g, added %g kcal/mol; dS = %G, added %g cal/mol;"
                       %(delta_h,tdh,delta_s,tds))
    if sequ.endswith('A'):
        tdh = nn_table['init_5T/A'][d_h]
        tds = nn_table['init_5T/A'][d_s]
        delta_h += tdh
        delta_s += tds
        if tdh != 0 or tds != 0:
            if ((disp_flag > 2)):
                print ("Penalty if 3' end is T dH = %g, added %g kcal/mol; dS = %G, added %g cal/mol;"
                       %(delta_h,tdh,delta_s,tds))
    if disp_flag > 2:
        checkdg = 0
        for i in range(0,nn-1):
            checkdg += dg[i]
        print ("\t1- dG = %g checkdg = %g" % (cor*(delta_h - ((delta_s) *(T_dG+273.15))/1000),checkdg))
    #print("after genral initiation   dg[%d] = %g kBT" % (nn-2,dg[nn-2]))
    # Type: Different values for G/C or A/T terminal basepairs
    ends = sequ[0] + sequ[-1]
    AT = ends.count('A') + ends.count('T')
    GC = ends.count('G') + ends.count('C')
    # in fact init_A/T = Init/2 + TerminalPenality and Init/2 = init_G/C (half the initiation penality)
    tdh = nn_table['init_A/T'][d_h]
    tds = nn_table['init_A/T'][d_s]
    delta_h += tdh * AT
    delta_s += tds * AT

    tmpdh = nn_table['init_A/T'][d_h] - nn_table['init_G/C'][d_h]
    tmpds = nn_table['init_A/T'][d_s] - nn_table['init_G/C'][d_s]
    dg_AT_penality = -cor*(tmpdh - (tmpds* TG/1000))

    if tdh * AT != 0 or tds * AT != 0:
        if ((disp_flag > 2)):
            print ("%d A/T terminal basepaairs dH += %g kcal/mol, dS += %g cal/mol, dH = %g; ds = %g;"
                   %(AT,tdh*AT,tds*AT,delta_h,delta_s))
    if (sequ[0] == 'A' or sequ[0] == 'T'):
        dg[0] += dg_AT_penality
        if ((disp_flag > 2)):
            print ("A/T terminal bp in 0     dH += %5.3f  kcal/mol dS += %5.3f  cal/mol dG += %5.3f  kcal/mol;\ndg[0] = %g kBT" %
                   (tmpdh, tmpds, dg_AT_penality,dg[0]))
    if (sequ[-1] == 'A' or sequ[-1] == 'T'):
        dg[nn-2] += dg_AT_penality
        if ((disp_flag > 2)):
            print ("A/T terminal bp in -1    dH += %5.3f kcal/mol dS += %5.3f cal/mol dG += %5.3f kcal/mol;\ndg[%d] = %g kBT" %
                   (tmpdh, tmpds, dg_AT_penality,nn-2,dg[nn-2]))
    #print("after A/T terminal bp in -1   dg[%d] = %g kBT" % (nn-2,dg[nn-2]))
    tdh = nn_table['init_G/C'][d_h]
    tds = nn_table['init_G/C'][d_s]
    delta_h += tdh * GC
    delta_s += tds * GC
    if tdh * GC != 0 or tds * GC != 0:
        if ((disp_flag > 2)):
            print ("%d G/C terminal basepairs  dH += %g kcal/mol, dS += %G cal/mol, dH = %g; ds = %g;"
                   %(GC,tdh*GC,tds*GC,delta_h,delta_s))
    if (sequ[0] == 'C' or sequ[0] == 'G'):
        if ((disp_flag > 2)):
            print ("G/C terminal bp in 0  no correction")
    if (sequ[-1] == 'C' or sequ[-1] == 'G'):
        if ((disp_flag > 2)):
            print ("G/C terminal bp in -1 no correction")


    if disp_flag > 2:
        tmpdh = nn_table['init_G/C'][d_h]
        tmpds = nn_table['init_G/C'][d_s]
        dgtt = -cor*(tmpdh - (tmpds* TG/1000))
        checkdg = 0
        for i in range(0,nn-1):
            checkdg += dg[i]
        print ("\t2- dG = %g checkdg = %g checkdg + 2*dgt = %g"
               % (cor*(delta_h - ((delta_s) *(T_dG+273.15))/1000),checkdg,checkdg+2*dgtt))

    # Finally, the 'zipping'
    for basenumber in range(len(tmp_seq) - 1):
        neighbors = tmp_seq[basenumber:basenumber + 2] + '/' + \
            tmp_cseq[basenumber:basenumber + 2]
        if neighbors in imm_table:
            tmpdh = imm_table[neighbors][d_h]
            tmpds = imm_table[neighbors][d_s]
            delta_h += tmpdh
            delta_s += tmpds
            dg[off_bp+basenumber] += -cor*(tmpdh - (tmpds * TG/1000))
            if ((disp_flag > 2)):
                print ("Zipping-%d of %s in mismatches\n\tdH += %5.2f kcal/mol, dS += %5.2f cal/mol, dH = %5.2f; ds = %5.2f;\n\tdg[%d] = %g kBT added %g kBT"
                       %(basenumber,neighbors,tmpdh,tmpds,delta_h,delta_s,off_bp+basenumber,dg[off_bp+basenumber],cor*(tmpdh - (tmpds * TG/1000))))
        elif neighbors[::-1] in imm_table:
            tmpdh = imm_table[neighbors[::-1]][d_h]
            tmpds = imm_table[neighbors[::-1]][d_s]
            delta_h += tmpdh
            delta_s += tmpds
            dg[off_bp+basenumber] += -cor*(tmpdh - (tmpds * TG/1000))
            if ((disp_flag > 2)):
                print ("Zipping-%d of %s in mismatches (%s) dH += %5.2f kcal/mol, dS += %5.2f cal/mol, dH = %5.2f; ds = %5.2f;\n\tdg[%d] = %g kBT  added %g kBT"
                       %(basenumber,neighbors[::-1],neighbors,tmpdh,tmpds,delta_h,delta_s,off_bp+basenumber,dg[off_bp+basenumber],cor*(tmpdh - (tmpds * TG/1000))))
        elif neighbors in nn_table:
            tmpdh = nn_table[neighbors][d_h]
            tmpds = nn_table[neighbors][d_s]
            delta_h += tmpdh
            delta_s += tmpds
            dg[off_bp+basenumber] += -cor*(tmpdh - (tmpds * TG/1000))
            if ((disp_flag > 2)):
                print ("Zipping-%d of %s\n\tdH += %5.2f kcal/mol, dS += %5.2f cal/mol, dH = %5.2f; ds = %5.2f;\n\tdg[%d] = %g kBT  added %g kBT"
                       %(basenumber,neighbors,tmpdh,tmpds,delta_h,delta_s,off_bp+basenumber,dg[off_bp+basenumber],cor*(tmpdh - (tmpds * TG/1000))))
        elif neighbors[::-1] in nn_table:
            tmpdh = nn_table[neighbors[::-1]][d_h]
            tmpds = nn_table[neighbors[::-1]][d_s]
            delta_h += tmpdh
            delta_s += tmpds
            dg[off_bp+basenumber] += -cor*(tmpdh - (tmpds * TG/1000))
            if ((disp_flag > 2)):
                print ("Zipping-%d of %s (%s) dH += %5.2f kcal/mol, dS += %5.2f cal/mol, dH = %5.2f; ds = %5.2f;\n\tdg[%d] = %g kBT added %g kBT"
                       %(basenumber,neighbors[::-1],neighbors,tmpdh,tmpds,delta_h,delta_s,off_bp+basenumber,dg[off_bp+basenumber],cor*(tmpdh - (tmpds * TG/1000))))
        else:
            # We haven't found the key...
            mt._key_error(neighbors, strict)

    #print("after zipping bp in -1   dg[%d] = %g kBT" % (nn-2,dg[nn-2]))
    # we assume that the hairpin has no mismatches and special ends!
    for basenumber in range(len(seq) - 1):
        hp_neighbors = h_seq[basenumber:basenumber + 2] + '/' + \
            c_seq[basenumber:basenumber + 2]
        if hp_neighbors in nn_table:
            tmpdh = nn_table[hp_neighbors][d_h]
            tmpds = nn_table[hp_neighbors][d_s]
            dgh[basenumber] += -cor*(tmpdh - (tmpds * TG/1000))
        elif hp_neighbors[::-1] in nn_table:
            tmpdh = nn_table[hp_neighbors[::-1]][d_h]
            tmpds = nn_table[hp_neighbors[::-1]][d_s]
            dgh[basenumber] += -cor*(tmpdh - (tmpds * TG/1000))
        else:
            # We haven't found the key...
            print("Cannot find %s in nn_table" % hp_neighbors)
            #mt._key_error(hp_neighbors, strict)



    if disp_flag > 2:
        tmpdh = nn_table['init_G/C'][d_h]
        tmpds = nn_table['init_G/C'][d_s]
        dgtt = -cor*(tmpdh - (tmpds* TG/1000))
        checkdg = 0
        for i in range(0,nn-1):
            checkdg += dg[i]
        print ("\t3- dG = %g checkdg = %g checkdg + 2*dgt = %g dgt = %g"
               % (cor*(delta_h - ((delta_s) *(T_dG+273.15))/1000),checkdg,checkdg+2*dgtt,dgtt))



    k = (dnac1 - (dnac2 / 2.0)) * 1e-9
    if selfcomp:
        k = dnac1 * 1e-9
        tdh = nn_table['sym'][d_h]
        tds = nn_table['sym'][d_s]
        delta_h += tdh
        delta_s += tds
        if tdh != 0 or tds != 0:
            if ((disp_flag > 2)):
                print ("Self complementary correction dH = %g, added %g kcal/mol; dS = %G, added %g cal/mol;"
                       %(delta_h,tdh,delta_s,tds))
    #R = 1.987  # universal gas constant in Cal/degrees C*Mol

    # We compute salt correction for the hybridized oligo that is with reduced charge near dsDNA
    corrh = salt_correction2(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs,
                            method=saltcorr, seq=seq)
    corr = dsCharge * corrh


    melting_temp = (1000 * delta_h) / (delta_s + (R * (math.log(k)))) - 273.15
    if saltcorr in (1, 2, 3, 4):
        melting_temp += corr
        delta_sc = ((1000 * delta_h) / (melting_temp + 273.15)) - (R * (math.log(k)))
        cordsh = delta_sc - delta_s
        cords = dsCharge * cordsh
    if saltcorr == 5 or saltcorr == 8:
        cords = corr
        cordsh = corrh
        melting_temp = (1000 * delta_h) / (delta_s + corr + (R * (math.log(k)))) - 273.15
        #print("M5 Tm = %g" % (melting_temp))
    if saltcorr in (6, 7):
        # Tm = 1/(1/Tm + corr)
        melting_temp = (1 / (1 / (melting_temp + 273.15) + corrh) - 273.15)
        delta_sc = ((1000 * delta_h) / (melting_temp + 273.15)) - (R * (math.log(k)))
        cordsh = delta_sc - delta_s
        cords = dsCharge * cordsh
        #print("M6 Tm = %g" % (melting_temp))
    for basenumber in range(len(seq) - 1):
        dg[basenumber] += cor*(cords/(nn-1))* TG/1000
        dgh[basenumber] += cor*(cords/(nn-1))* TG/1000 # was cor*(cordsh/(nn-1))* TG/1000
        if len(energy_mod) >= (len(seq) - 1):
            dg[basenumber] += energy_mod[basenumber]
            if ((disp_flag > 2)):
                  print("dg[%d] modified by %8.5f" % (basenumber,energy_mod[basenumber])) 
        if ((disp_flag > 2)):
            print("dg[%d] = %8.5f kBT; dgh[%d] = %8.5f kBT += %8.5f kBt" %
                  (basenumber,dg[basenumber],basenumber,dgh[basenumber],cor*(cords/(nn-1))* TG/1000))



    # we compute salt correction for normal oligo to compute Kd and encercling energy
    if saltcorr:
        corr = salt_correction2(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs,
                               method=saltcorr, seq=seq)
        cordsbp = corr/(nn-1)


    #Delta_G = RT ln(Kd)
    # exp(Delta_G/RT) = Kd
    # exp( (delta_h - ((delta_s+corr) *(T_dG+273.15))/1000)RT) = Kd
    #
    #if saltcorr == 5:
    #    delta_g = delta_h - ((delta_s+corr) *(T_dG+273.15))/1000;
    #    melting_temp = (1000 * delta_h) / (delta_s + corr + (R * (math.log(k)))) - 273.15
    #else:
    #    delta_g = delta_h - ((delta_s) *(T_dG+273.15))/1000;
    #    melting_temp = (1000 * delta_h) / (delta_s + (R * (math.log(k)))) - 273.15

    melting_temp = (1000 * delta_h) / (delta_s + (R * (math.log(k)))) - 273.15
    if saltcorr in (1, 2, 3, 4):
        melting_temp += corr
        delta_sc = ((1000 * delta_h) / (melting_temp + 273.15)) - (R * (math.log(k)))
    if saltcorr == 5 or saltcorr == 8:
        #tmp = (delta_h * 1000)/(T_dG+273.15);
        #tmp -= (delta_s + corr);
        #tmp /= R;
        #tmp = numpy.exp(tmp);
        #Kd = tmp * 1000000
        delta_sc = delta_s + corr
        melting_temp = (1000 * delta_h) / (delta_s + corr + (R * (math.log(k)))) - 273.15
    if saltcorr in (6, 7):
        # Tm = 1/(1/Tm + corr)
        delta_g = delta_h - ((delta_s) *(T_dG+273.15))/1000;
        melting_temp = (1 / (1 / (melting_temp + 273.15) + corr) - 273.15)
        delta_sc = ((1000 * delta_h) / (melting_temp + 273.15)) - (R * (math.log(k)))
        #tmp = (delta_h * 1000)/(T_dG+273.15);
        #tmp -= (delta_sc);
        #tmp /= R;
        #tmp = numpy.exp(tmp);

    delta_g = delta_h - ((delta_sc) *(T_dG+273.15))/1000;
    delta_g37 = delta_h - ((delta_sc) *(37+273.15))/1000;
    cordsbp = (delta_sc - delta_s)/(nn-1)

    roo = numpy.zeros(nn)
    roo2 = numpy.zeros(nn)
    roh = numpy.zeros(nn)
    tmpdh = nn_table['init_G/C'][d_h]
    tmpds = nn_table['init_G/C'][d_s]
    # In fact init_G/C equals half the initiation energy 
    dg_half_init = -cor*(tmpdh - (tmpds* TG/1000))

    checkdg = 0

    for i in range(0,nn-1):
        neighbors = tmp_seq[i:i+2] + '/' + tmp_cseq[i:i+2]
        roo[i]=numpy.exp(-dg[i])
        roh[i]=numpy.exp(-dgh[i])
        # we decrease every energy by E_init/2, thus the energy is reduced for the last two bases hybridized
        roo2[i]=numpy.exp(-dg[i]-dg_half_init)
        checkdg += dg[i]
        if disp_flag > 2:
            print ("For %s: dg = %8.5f kBT dgh = %8.5f kBT dg2 = %8.5f kBT" % (neighbors,dg[i],dgh[i],dg[i]+dg_half_init))
    checkdg += 2*dg_half_init

    if disp_flag > 2:
        print ("\t4- dG = %g checkdg (+2 half_init) = %g dg_half_init = %g " % (cor*delta_g,checkdg, dg_half_init))

    '''
    tmpdh = nn_table['init_A/T'][d_h] - nn_table['init_G/C'][d_h]
    tmpds = nn_table['init_A/T'][d_s] - nn_table['init_G/C'][d_s]
    dg_AT_penality = -cor*(tmpdh - (tmpds* TG/1000))
    if seq[0] == 'A' or  seq[0] == 'T' or seq[0] == 'a' or  seq[0] == 't':
        roo[0]=numpy.exp(-dg[0] - dg_AT_penality)
        roo2[0]=numpy.exp(-dg[0] - dg_AT_penality - dg_half_init)
        checkdg += dg_AT_penality
        if ((disp_flag > 2)):
            print ("dg[0] = %8.5f kBT dg2 = %8.5f kBT"
                   % (dg[0] + dg_AT_penality ,dg[0] + dg_AT_penality + dg_half_init))

    if seq[-1] == 'A' or  seq[1] == 'T' or seq[-1] == 'a' or  seq[-1] == 't':
        roo[nn-2]=numpy.exp(-dg[nn-2] + dg_AT_penality)
        roo2[nn-2]=numpy.exp(-dg[nn-2] + dg_AT_penality - dg_half_init)
        checkdg -= dg_AT_penality
        if disp_flag > 2:
            print ("dg[%d] = %8.5f kBT dg2 = %8.5f kBT"
                   % (nn-2,dg[nn-2] - dg_AT_penality,dg_half_init - dg_AT_penality))
    '''
    if disp_flag > 2:
        print ("\t5- dG = %g checkdg (+2 half_init) = %g dg_half_init = %g " % (cor*delta_g,checkdg, dg_half_init))


# roo[0] opening of the first 5' base of the oligo towards the fork

    if disp_flag > 2 :
        for i in range(0,nn-1):
            print("Roh[%d] = %8.5f; roo[%d] = %8.5f; roo2[%d] = %8.5f"
                  % (i,roh[i],i,roo[i],i,roo2[i]));

#elasticity models for the closing rates for the ds and ss DNA
    #d = 0.542 # nm, 0.567 #0.56*1.12 new 0.605
    #bp = 2.14# nm,1.828 # 1.78 # 1.5*1.125 new 1.828
    #S = 216# pN

    d = dss
    bp = bpss
    S = Sss

    ft = 4.1*(T0+T_dG)/(T0+25);

    u = bp*F/ft;
    if (u == 0) :
        gss = 0
    else:
        gss = d/bp*(numpy.log(numpy.sinh(u)/u)+ 0.5*F**2./ft/S);
    gds=(F/ft - numpy.sqrt(F/ft/50) + 0.5*F**2./ft/1230)*0.34;

#closing rates:only depend on the force
    rch = numpy.exp(-2*gss)
    rco = numpy.exp(-gss+gds)
# encircling rate for one bounded bp of the oligo
    tmpdg = -2*gss -gds
    if ((disp_flag > 1) and (loop == 1)):
        print("encercling energy %g" % (tmpdg));
    # must depend on salt corr
    # was tmpdg += 3*(0.368*numpy.log(0.001*Na+3.2*numpy.sqrt(0.001*Mg))) * (T_dG+T0)/1000;
    tmpdg += encerclingBpSalt * cordsbp * (T_dG+T0)/1000;
    if ((disp_flag > 1) and (loop == 1)):
        print("encercling energy cor %g" % (tmpdg));
    ren = encerclingCoeff*numpy.exp(tmpdg)  # was 0.1
    if encercling == 0:
        ren = 0
    if (disp_flag > 1):
        print("final encercling energy rate %g encercling %d loop %d" % (ren,encercling,loop));        
    if (disp_flag > 1):
        if (loop == 0):
            print("Hairpin closing rate = %g oligo closing rate = %g (depend on force only F = %g)" % (rch,rco,F));
        else :
            print("oligo closing rate = %g (depends on force only F = %g)" % (rco,F));

# l is an index over all possible configurations
#ind[i, j]  is the index of the configuration :
# i = position of the opening front on the 5' end, [0, nn-1]
# j = position of the opening front on the 3' end, nn-1 nothing open

    n = nn-1;
    l=0
# loop version
    if (loop == 1) :
        ind=numpy.zeros((nn, nn),dtype=numpy.int32)
        for i in range(0,nn):
            for j in range (i+1,nn):
                ind[i,j]=l
                l=l+1
        n = l;
        if ((disp_flag > 3) and (nn < 10)):
            print(" \nindex over all possible configurations ind[Oligo 5' end position, oligo 3'end position]");
            for i in range(0,nn):
                for j in range(i+1,nn):
                    print("ind[%d,%d] = %g" % (i,j,ind[i,j]));
# transition matrix
        mat=numpy.zeros((l,l))
        mat = fill_loop_mat(mat,ind,l,nn,rco,roo,roo2,ren,disp_flag,LAST_2bp)
        if (disp_flag >= 3):
            print("before ------------------------------------- \n")
            print("-m 1\n-a \n-src \"roh data over %d points, rco = %g, rch = %g Na = %g K = %g Tris = %g Mg = %g dNTPs = %g saltcor %d T = %g dsCharge = %g F = %g rate = %g \""% (nn,rco,rch,Na, K, Tris, Mg, dNTPs, saltcorr, T_dG, dsCharge, F, rate));
            for i in range(0,nn-1):
                print("%g\t" % (roh[i]), end='', flush=True);
            print("\n-m 1\n-a \n-src \"roo data over %d points, rco = %g, rch = %g\" \n"% (nn,rco,rch));
            for i in range(0,nn-1):
                print("%g\t" % (roo[i]), end='', flush=True);
            print("\n-m 1\n-a \n-src \"roo2 data over %d points, rco = %g, rch = %g\" \n"% (nn,rco,rch));
            for i in range(0,nn-1):
                print("%g\t" % (roo2[i]), end='', flush=True);                
            print("\nafter------------------------------------ \n")
        
    elif (loop == 0) :
        ind=numpy.zeros((nn, nn, nn),dtype=numpy.int32)
        for i in range(0,nn):
            for j in range (0,i+1):
                for k in range (i+1,nn):
                    ind[i,j,k]=l
                    l=l+1

        n = l;
        if ((disp_flag > 3) and (nn < 8)):
            disp3_i = 0
            mess3 = ""
            d1 = [-1,-1,-1]
            d2 = [-1,-1,-1]
            d3 = [-1,-1,-1]
            print(" \nWe describe below index over all the possible configurations: ind[5'oligo position, fork, 3'oligo position]");
            for i in range(0,nn):
                for j in range(0,i+1):
                    for k in range (i+1,nn):
                        mess3 = mess3 +("      ind[%d,%d,%d] = %g      " % (i,j,k,ind[i,j,k]));
                        d1[disp3_i] = i
                        d2[disp3_i] = j
                        d3[disp3_i] = k
                        disp3_i = disp3_i + 1
                        if disp3_i > 2 :
                            print(mess3)
                            disp_three_fork(nn-1,d1[0],d2[0],d3[0],d1[1],d2[1],d3[1],d1[2],d2[2],d3[2])
                            disp3_i = 0
                            mess3 = ""
                            d1 = [-1,-1,-1]
                            d2 = [-1,-1,-1]
                            d3 = [-1,-1,-1]
            if disp3_i > 0 :
                print(mess3)
                disp_three_fork(nn-1,d1[0],d2[0],d3[0],d1[1],d2[1],d3[1],d1[2],d2[2],d3[2])

# transition matrix
        mat=numpy.zeros((l,l))
        if (disp_flag >= 3):
            print("before ------------------------------------- \n")
            print("-m 1\n-a \n-src \"roh data over %d points, rco = %g, rch = %g Na = %g K = %g Tris = %g Mg = %g dNTPs = %g saltcor %d T = %g dsCharge = %g F = %g rate = %g \""% (nn,rco,rch,Na, K, Tris, Mg, dNTPs, saltcorr, T_dG, dsCharge, F, rate));
            for i in range(0,nn-1):
                print("%g\t" % (roh[i]), end='', flush=True);
            print("\n-m 1\n-a \n-src \"roo data over %d points, rco = %g, rch = %g\" \n"% (nn,rco,rch));
            for i in range(0,nn-1):
                print("%g\t" % (roo[i]), end='', flush=True);
            print("\n-m 1\n-a \n-src \"roo2 data over %d points, rco = %g, rch = %g\" \n"% (nn,rco,rch));
            for i in range(0,nn-1):
                print("%g\t" % (roo2[i]), end='', flush=True);                
            print("\nafter------------------------------------ \n")
        mat = fill_fork_mat(mat,ind,l,nn,rco,rch,roo,roo2,roh,disp_flag,LAST_2bp)


#holding state: all possible configurations
    un=numpy.ones(n)

#initial state: nn open bases for the hairpin 0 open bases for the oligo
    p0=numpy.zeros(n)
    if (loop == 1) :
        for i in range(0,int(l)):
            if(i == ind[0,nn-1]):
                p0[i]=1
        if ((disp_flag > 1)):
            print ("Oligo on template apex at F= %g" % (F))

    elif (loop == 0) :
        for i in range(0,int(l)):
            if(i == ind[0,0,nn-1]):
                p0[i]=1
        if ((disp_flag > 1)):
            print ("Oligo displeced by fork at F= %g" % (F))


    if (disp_flag > 3) and ( nn < 8):
        print("\nP0 = \n");
        print(p0);
        print("\nMatrice\n");
        print(mat);
    p0T= numpy.matrix(p0).T
#inverse of the transition matrix

    if (disp_flag > 3) and ( nn < 8):
        print("\nTransposed Matrix\n");
        print(p0T);
    matinv=numpy.linalg.inv(mat)
    if (disp_flag > 3) and ( nn < 8):
        print("\nInvert Matrix\n");
        print(matinv);
# between the initil and final state
    ra=rate*10**(-6.); #1.33*10**(-6.);
    trep=-ra*numpy.matrix(un*(numpy.matrix(matinv*p0T)))

    dgf = cor*delta_g+((nn-1)*(gss-gds))
    tmp = numpy.exp(dgf);
    Kd = tmp * 1000000
    koff = 1/trep
    kon = koff/Kd
    template = ''
    if c_seq_defined > 0:
        template = " on template 3'-%s-5'" % (c_seq)
    if ((disp_flag > 0)):
        #print ("For salt correction %d, cordsbp %g" % (saltcorr,cordsbp))
        print ("\n-m 1\n-src \"For %s oligo 5'-%s-3' %d bases%s,%s T=%gC F=%gpN; Na+=%gmM; Mg++=%gmM; tris=%gmM; dsCharge=%g; rate=%g"%(oligotype,seq,len(seq),template,mod_erg,T_dG,F,Na,Mg,Tris,dsCharge,ra))
        print ("dH = %g Kcal/mol; dS0 = %g; dS = %g (cal/mol);\ndG(%gC) = %g (kcal/mol) = %g kBT; dG(37C) = %g (kcal/mol) Kd = %g micoM" %(delta_h,delta_s,delta_sc,T_dG,delta_g,cor*delta_g,delta_g37,Kd))
        if selfcomp == True :
            print ("Tm(F=0) = %g T = %g; Na+ = %g; Mg2+ = %g; dnac1 = %g nM; dnac2 = %g nM; Salt Correction %d self-complenentary" % (melting_temp,T_dG,Na,Mg,dnac1,dnac2,saltcorr))
        else :
            print ("Tm(F=0) = %g T = %g; Na+ = %g; Mg2+ = %g; dnac1 = %g nM; dnac2 = %g nM; Salt Correction %d non self-complenentary" % (melting_temp,T_dG,Na,Mg,dnac1,dnac2,saltcorr))
            print ("Kon = %g x 10^6; Koff = %g s-1 Kd = %g Toff %g S; dG(F=%gpN) = %g kBT" % (kon[0,0],koff[0,0],Kd,trep[0,0],F,dgf))
        print("ssDNA elasticity d = %g nm; bp = %g nm; S = %g pN" % (dss,bpss,Sss))
        print("at F = %g pN rezipping energy/bp = %g kBT Dg_{avg}/bp = %g kBt" % (F,2*gss,cor*delta_g/(len(seq)-1)))
        if (loop == 1) and (encercling == 1) :
            print ("Oligo in the apex of an hairpin");
            print ("For encercling, general rate coefficient = %g, entropy salt correction extend in nts = %g\"" % (encerclingCoeff,encerclingBpSalt))
        if (loop == 1) and (encercling == 0) :
            print ("Simple Oligo hybridization\"");
        if (loop == 0):
            print ("Oligo blocking a fork rezipping\"");
        #print("Reduce energy by Initiation on last two bases bounded = %d\n"
        #      "Reduce energy on 5' by AT penalty = %d\n"
        #      "Reduce energy on 3' by AT penalty = %d\n" % (reduce_by_IN_on_last_2_bases
        #                                                    ,use_AT_penality_5
        #                                                    ,use_AT_penality_3))

    #tmp = 0.000001*oligoConc*2; # we are at 1 microM
    #tmp = R * numpy.log(tmp/4); # 4 for non self complementary
    #tmp += S98;
    #tmp = 1000/tmp;
    #Tm98 = (H98 * tmp) - T0;


    return trep,koff,kon,dgf, melting_temp



#mystring = 'GGGTTG'  #'ACAAGTCCT'
#mystring = 'AAGTCCT'
#myseq = Seq.Seq(mystring)
#print('%0.2f' % mt.Tm_Wallace(mystring))
#    84.00
#print('%0.2f' % mt.Tm_Wallace(myseq))
#    84.00
#print('%0.2f' % mt.Tm_GC(myseq))
#    58.73
#print('%s Tm=%0.2f' % (mystring,Tm_NN2(myseq,Na=125,Mg=0,dnac1=200, dnac2=200,T_dG=22, saltcorr=5, encerclingCoeff=0, f=0, encercling=0,koffNa=100)))
#print('%s Tm=%0.2f' % (mystring,Tm_NN2(myseq,Na=125,Mg=0,dnac1=200, dnac2=200,T_dG=22, saltcorr=6, encerclingCoeff=0, f=0, encercling=0,koffNa=100)))
#    60.32 ,c_seq="TGTTCAGGA" c_seq="TGTTCAGGA",


#Ha oligos "TGTTCAGGA"
#first 'ACAAGTCCT' Kon = 0.316937 x 10^6; Koff = 0.111438 Kd = 0.351611
#mis 1 'TCAAGTCCT' Kon = 0.244412 x 10^6; Koff = 0.256973 Kd = 1.05139
#mis 2 'AGAAGTCCT' Kon = 0.0324713 x 10^6; Koff = 7.16304 Kd = 220.596
#mis 3 'ACTAGTCCT' Kon = 0.142446 x 10^6; Koff = 16.91 Kd = 118.712
#mis 4 'ACATGTCCT' Kon = 0.193367 x 10^6; Koff = 12.8303 Kd = 66.3521
#mis 5 'ACAACTCCT' Kon = 0.0351204 x 10^6; Koff = 305.375 Kd = 8695.1
#mis 6 'ACAAGACCT' Kon = 0.193356 x 10^6; Koff = 21.8126 Kd = 112.81
#mis 7 'ACAAGTGCT' Kon = 0.195698 x 10^6; Koff = 10.2569 Kd = 52.4119
#mis 8 'ACAAGTCGT' Kon = 0.0903952 x 10^6; Koff = 5.49138 Kd = 60.7486
#mis 9 'ACAAGTCCA' Kon = 0.220663 x 10^6; Koff = 0.166578 Kd = 0.754898






def usage():
    print ("This program computes the displacement time of an oligonuleotide in the apex of a hairpin", file=sys.stderr)
    print ("The algorithm uses (Initiation Energy) for the last two bases and AT penalty at the oligo ends", file=sys.stderr)
    print ("It uses the elastic FJC model for ssDNA elasricity", file=sys.stderr)
    print ("pushed by a DNA fork\nYou must specify the oligonucleotide sequence,\nyou may specify the force and the temperature", file=sys.stderr)
    print ("The program accepts GNU type options:", file=sys.stderr)
    print ("\nExample 1: Tdisp -O=ACGTTCGA", file=sys.stderr)
    print ("displays Tdisp for the oligo with the sequence ACGTTCGA \nat force = 8.5 pN and T = 25 C", file=sys.stderr)
    print ("the oligonucleotide sequence may be specified using: -O=... or --oligo=...", file=sys.stderr)
    print ("\nExample 2: Tdisp-Fork-Blocking -O=ACGTTCGA -F=10.2 -T=30.5", file=sys.stderr)
    print ("will display Tdisp for the oligo having the sequence ACGTTCGA \nat force = 10.2 pN and T = 30.5 C", file=sys.stderr)
    print ("the force may be specified using: -F=x or --force=x ", file=sys.stderr)
    print ("the temperature may be specified using: -T=t or --temperature=x ", file=sys.stderr)
    print ("You can scan force using -R=fmin:fmax:fstep or --forceRange=... ", file=sys.stderr)
    print ("The order of options is irrelevant, force and temperature are optional", file=sys.stderr)
    print ("-O or --oligo            => define oligo : -O=ACGTA", file=sys.stderr)
    print ("      --oligotype        => may be DNA, RNA or LNA", file=sys.stderr)
    print ("      --template         => define template --template=TGCAT", file=sys.stderr)
    print ("      --Mg               => specify [mg2+] in mM (default 0 mM)", file=sys.stderr)
    print ("      --Na               => specify [mg2+] in mM (default 150 mM)", file=sys.stderr)
    print ("      --tris             => specify [Tris] in mM (default 20 mM)", file=sys.stderr)
    print ("      --oligoConc        => specify [oligo] in nM (default 25 nM)", file=sys.stderr)
    print ("-T or --temperature      => define temperature -T=25", file=sys.stderr)
    print ("-F    --force            => specify Force in pN -F=8 ", file=sys.stderr)
    print ("       --mode            => either : oligo (default), fork or apex", file=sys.stderr)
    print ("-R    --forceRange       => iterate in force -F=2:12:1 (from 2 pN to 12 in step of 1", file=sys.stderr)
    print ("      --tempRange        => iterate in temperature --tempRange=20:32:1 (from 20 degre to 32 in step of 1", file=sys.stderr)
    print ("-v                       => verbose", file=sys.stderr)
    print ("-h or --help             => help", file=sys.stderr)
    print ("      --output Toff/koff/kon/dGF ", file=sys.stderr)
    print ("      --rate             => intrinsic ration rate (in micro seconds)", file=sys.stderr)
    print ("      --dsCharge         => modifier of dsDNA salt correction", file=sys.stderr)
    print ("-D    --debug            => Debug level", file=sys.stderr)
    print ("      --grout            => specify a gr output for graphics ", file=sys.stderr)
    print ("      --salt             => use salt correction --sal=5 or --salt=6 (default)", file=sys.stderr)
    print ("      --loop             => specify that this is an hairpin substrate", file=sys.stderr)
    print ("      --kon              => give kon as output (instead of Toff)", file=sys.stderr)
    print ("      --midInitiation    => reduce the energy of the last 2 bases attached by 1/2 of initiation energy", file=sys.stderr)
    print ("      --AT_penality_5    => use AT penalty on 5'", file=sys.stderr)
    print ("      --AT_penality_3    => use AT penalty on 3'", file=sys.stderr)
    print ("      --ssDNA_d          => specify the ssDNA monomer length (0.54nm)", file=sys.stderr)
    print ("      --ssDNA_bp         => specify the ssDNA persistence length (2.14nm)", file=sys.stderr)
    print ("      --ssDNA_S          => specify the ssDNA enthalpic term (216 pN)", file=sys.stderr)
    print ("      --modify           => do not use", file=sys.stderr)
    print ("-E    --encercling       => allow encercling -E=1 or not -E=0", file=sys.stderr)
    print ("      --encerclingCoeff  => define encercling coefficient ", file=sys.stderr)
    print ("      --encerclingBpSalt => specify the nb. of Bp in encercling", file=sys.stderr)




def main():
    seq = ""
    T = 25;
    f = 8.5;
    debug = 1;
    fmin = -1;
    fmax = -1;
    fstep = -1;
    Tmin = -1;
    Tmax = -1;
    Tstep = -1;

# default settings
    encercling_enable = 0
    reduce_by_IN_on_last_2_bases = 1
    #rate0 = 1 #4/1.4
    rate0 = 2.5 
    use_salt_correction = 6
    use_AT_penality_5 = 1
    use_AT_penality_3 = 1
    Na = 150
    Mg = 0
    tris = 20
    dsCharge = 0.5
    oligoConc = 25
    encerclingCoeff = 0.31
    encerclingBpSalt = 2
    #loop = 1
    loop = 0
    kon = 0
    oligotype ='DNA'
    modify = ""
    debug = 1
    mode = ''
    Output = 'Toff'
    Grout = 0
    ssDNA_d = 0.542
    ssDNA_bp = 2.14
    ssDNA_S = 216
    fork_tau_mul = 2.5

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:O:T:F:R:D:E:v",
                                   ["help", "oligo=","temperature=","force=","forceRange="
                                    ,"debug=","encercling=","salt=","rate=","midInitiation="
                                    ,"Na=","Mg=","AT_penality_5=","AT_penality_3=","loop="
                                    ,"oligoConc=","encerclingCoeff=","kon=","dsCharge="
                                    ,"tempRange=","oligotype=","verbose","modify="
                                    ,"template=","mode=","tris=","output=","encerclingBpSalt="
                                    ,"grout=","ssDNA_d=","ssDNA_bp=","ssDNA_S="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print (err, file=sys.stderr)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    seq = "ATCGAATCGA"
    C_seq = None
    verbose = False
    F = 8.5;
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-O"):
            seq = a[1:]
        elif o in ("--oligo"):
            seq = a[0:]
        elif o in ("--template"):
            C_seq = a[0:]
        elif o in ("--output"):
            Output = a[0:]
        elif o in ("--mode"):
            mode = a[0:]
            if (mode == 'oligo'):
                if (debug > 1):
                    print("mode oligo")
                encercling_enable = 0
                loop = 1
                if (f < 0):
                    f = 0
            elif (mode == 'fork'):
                if (debug > 1):
                    print("mode fork")
                loop = 0
                if (f < 0):
                    f = 8.5
                rate0 = rate0 * fork_tau_mul
            elif (mode == 'apex'):
                if (debug > 1):
                    print("mode apex")
                encercling_enable = 1
                if (f < 0):
                    f = 5
                loop = 1
        elif o in ("-T"):
            T = float(a[1:])
        elif o in ("--temperature"):
            T = float(a[0:])
        elif o in ("--dsCharge"):
            dsCharge = float(a[0:])
        elif o in ("-D"):
            debug = int(a[1:])
        elif o in ("--debug"):
            debug = int(a[0:])
        elif o in ("--grout"):
            Grout = int(a[0:])
        elif o in ("-E"):
            encercling_enable = int(a[1:])
        elif o in ("--encercling"):
            encercling_enable = int(a[0:])
        elif o in ("--encerclingBpSalt"):
            encerclingBpSalt = float(a[0:])
        elif o in ("--salt"):
            use_salt_correction = int(a[0:])
            if use_salt_correction == 6:
                #rate0 = rate0 * 0.8
                dsCharge = 0.5
                encerclingBpSalt = 2
                encerclingCoeff = 0.27
            elif use_salt_correction == 5 or use_salt_correction == 8:
                dsCharge = 0.5
                encerclingBpSalt = 2
                encerclingCoeff = 0.31
        elif o in ("--loop"):
            loop = int(a[0:])
        elif o in ("--oligotype"):
            oligotype = a[0:]
        elif o in ("--kon"):
            kon = int(a[0:])
        elif o in ("--midInitiation"):
            reduce_by_IN_on_last_2_bases = int(a[0:])
        elif o in ("--Mg"):
            Mg = float(a[0:])
        elif o in ("--Na"):
            Na = float(a[0:])
        elif o in ("--tris"):
            tris = float(a[0:])
        elif o in ("--oligoConc"):
            oligoConc = float(a[0:])
        elif o in ("--AT_penality_5"):
            use_AT_penality_5 = int(a[0:])
        elif o in ("--AT_penality_3"):
            use_AT_penality_3 = int(a[0:])
        elif o in ("-F"):
            f = float(a[1:])
        elif o in ("--force"):
            f = float(a[0:])
        elif o in ("--ssDNA_d"):
            ssDNA_d = float(a[0:])
        elif o in ("--ssDNA_bp"):
            ssDNA_bp = float(a[0:])
        elif o in ("--ssDNA_S"):
            ssDNA_S = float(a[0:])
        elif o in ("--modify"):
            modify = a[0:]
        elif o in ("--encerclingCoeff"):
            encerclingCoeff = float(a[0:])
            encercling_enable = 1
        elif o in ("--encerclingBpSalt"):
            encerclingBpSalt = float(a[0:])
            encercling_enable = 1
        elif o in ("--rate"):
            rate0 = float(a[0:])
        elif o in ("-R"):
            (fmin, fmax, fstep) = [t(s) for t,s in zip((float,float,float),a[1:].split(":"))]
        elif o in ("--forceRange"):
            (fmin, fmax, fstep) = [t(s) for t,s in zip((float,float,float),a[0:].split(":"))]
        elif o in ("--tempRange"):
            (Tmin, Tmax, Tstep) = [t(s) for t,s in zip((float,float,float),a[0:].split(":"))]
        elif o in ("--verbose"):
            print ("Options")
            print ("Force = %g pN" % (f))
            print ("Temperature = %g degrees C" % (T))
            print ("Oligo %s type %s" % (seq,oligotype))
            print ("Oligo concentration to compute Tm %g micoMolar" % (oligoConc))
            print ("correction by [Na+] in solution %g mM" % (Na))
            print ("correction by [Mg2+] in solution %g mM" % (Mg))
            print ("salt correction near dsDNA is reduced by %g" % (dsCharge))
            print ("Geometry of the experiment:")
            if (loop > 0 and encercling > 0):
                geo = "Oligo blocking hairpin in the apex"
            elif (loop > 0 and encercling == 0):
                geo = "Oligo binding a complementary template"
            elif (loop == 0):
                geo = "Oligo blocking hairpin pushed by the refolding fork"
            print ("%s" % (geo))
            if (loop > 0):
                print ("For oligo in the apex of a Hairpin:")
                print ("Encercling coeff %g, Nb. of bases impacr by salt in encercling %g"
                       % (encerclingCoeff,encerclingBpSalt))
            print ("dG correction by initiation and AT penality:")
            print ("Reduce energy on the last two bases by initiation %d\n"
                   "AT penalty on 5'%d\nAT penalty on 3'%d"
                   % (reduce_by_IN_on_last_2_bases,use_AT_penality_5
                      ,use_AT_penality_3))
            print ("MT timing %g mico secondes" % (rate0))
            return
        else:
            print ("Wrong option %s!" % (o), file=sys.stderr)
            assert False, "unhandled option"
    if len(seq) < 1:
        print ("Wrong seqence %s" % (seq), file=sys.stderr)  # will print something like "option -a not recognized"
        sys.exit(2)


    # ...
    #print(options)
    #print ("Fmin %g fmax= %g fstep = %g" % (fmin,fmax,fstep))
    #if (fmin < 0 ) and (fmax < 0) and (fstep < 0) and (Tmin < 0) and (Tmax < 0) and (Tstep < 0):
        #trep1 = compute_toff(seq, T, f, encercling_enable, reduce_by_IN_on_last_2_bases, rate0, use_salt_correction,
        #      use_AT_penality_5, use_AT_penality_3, Na, Mg, koffNa, koffMg,
        #       koffSaltCoeff, oligoConc, encerclingCoeff, encerclingBpSalt,
        #       loop, kon, oligotype, modify, debug )
        #print("%s\t%g" % (seq,trep1));
    #el
    NN_table=mt.DNA_NN3
    NN_table_desc = "DNA NN : Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594"
    TMM_table=mt.DNA_TMM1
    TNN_table_desc = ("# Terminal mismatch table (DNA)\n"
                      "# SantaLucia & Peyret (2001) Patent Application WO 01/94611\n")
    IMM_table=mt.DNA_IMM1
    INN_table_desc = ("# Internal mismatch and inosine table (DNA)\n"
                      "# Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594\n"
                      "# Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444\n"
                      "# Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179\n"
                      "# Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701\n"
                      "# Peyret et al. (1999), Biochemistry 38: 3468-3477\n"
                      "# Watkins & SantaLucia (2005), Nucl Acids Res 33: 6258-6267\n")
    DE_table=mt.DNA_DE1
    DE_table_desc = ("# Dangling ends table (DNA)\n"
                     "# Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934")




    if (oligotype == 'RNA'):
        if (debug > 1):
            print("RNA oligo")
        NN_table = RNA_DNA_NN1 #mt.R_DNA_NN1
        NN_table_desc = ("# RNA/DNA\n"
                         "# Sugimoto et al. (1995), Biochemistry 34: 11211-11216\n")

        TMM_table = mt.DNA_TMM1
        IMM_table = mt.DNA_IMM1
        DE_table = mt.DNA_DE1,
        DE_table_desc = ("RNA # Dangling ends table (RNA)\n"
                         "# Turner & Mathews (2010), Nucl Acids Res 38: D280-D282\n")
    elif (oligotype == 'LNA'):
        if (debug > 1):
            print("LNA oligo")
        NN_table = LNA_DNA_NN1
        NN_table_desc = ("# Owczarzy R, You Y, Groth CL, and Tataurov AV (2011) Stability and mismatch discrimination of locked nucleic acid-DNA duplexes. Biochemistry, 50:9352-9367.")
        TMM_table = LNA_DNA_TMM1 #mt.DNA_TMM1
        IMM_table = LNA_DNA_IMM1
        INN_table_desc = ("# LNA/DNA\n"
                         "# Owczarzy R, You Y, Groth CL, and Tataurov AV (2011) Stability and mismatch discrimination of locked nucleic acid-DNA duplexes. Biochemistry, 50:9352-9367.")

        DE_table = mt.DNA_DE1,
    elif (oligotype == 'LNA2'):
        if (debug > 1):
            print("LNA2 oligo")
        NN_table = LNA_DNA_NN2
        NN_table_desc = ("# modified from Owczarzy R, You Y, Groth CL, and Tataurov AV (2011) Stability and mismatch discrimination of locked nucleic acid-DNA duplexes. Biochemistry, 50:9352-9367.")
        TMM_table = mt.DNA_TMM1
        IMM_table = LNA_DNA_IMM1
        INN_table_desc = ("# LNA/DNA\n"
                         "# Owczarzy R, You Y, Groth CL, and Tataurov AV (2011) Stability and mismatch discrimination of locked nucleic acid-DNA duplexes. Biochemistry, 50:9352-9367.")
        DE_table = mt.DNA_DE1,
    elif (oligotype == 'DNA'):
        if (debug > 1):
            print("DNA oligo")
        NN_table = mt.DNA_NN3
        TMM_table = mt.DNA_TMM1
        IMM_table = mt.DNA_IMM1
        DE_table = mt.DNA_DE1,
    elif (oligotype == 'DNA2'):
        if (debug > 1):
            print("DNA2 oligo")
        NN_table = DNA_DNA_NN4
        TMM_table = mt.DNA_TMM1
        IMM_table = mt.DNA_IMM1
        DE_table = mt.DNA_DE1,
    else :
        print("wrong oligo type %s" % oligotype)


    if ((debug > 2)):
        print("\nWe are going to compute the free energy profile using nearest neighbour algorithm relying on:\n "
        "https://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html \nwritten by  Sebastian Bassi and Markus Piotrowski.\n"
              "We compute both the total free energy dG and the profile dg[].\n"
              "For DNA by default the energies of dinucleotides are taken from Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594\n")
        if (use_salt_correction == 5) :
            print("Salt correction is done using method 5 of MetltTemp.py using Correction for deltaS: 0.368 x (N-1) x ln[Na+] (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465)")
        if (use_salt_correction == 8) :
            print("Salt correction is done using method 5 of MetltTemp.py using Correction for deltaS: 0.368 x (N) x ln[Na+] (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465) to check other prog")
        if (use_salt_correction == 6) :
            print("Salt correction is done using method 6 of MetltTemp.py using (4.29(%GC)-3.95)x1e-5 x ln[Na+] + 9.40e-6 x ln[Na+]^2 Owczarzy et al. (2004), Biochemistry 43: 3537-3554\n)")
    if (debug > 2) :
        print(NN_table_desc)
        if (debug > 6) :
            print(NN_table)
        if (C_seq != None) :
            print(INN_table_desc)
            if (debug > 6) :
                print(INN_table)
            print(TNN_table_desc)
            if (debug > 6) :
                print(TNN_table)
            print(DE_table_desc)
            if (debug > 6) :
                print(DE_table)


    if (fmin < 0 ) and (fmax < 0) and (fstep < 0) and (Tmin < 0) and (Tmax < 0) and (Tstep < 0):
        trep1 = Tm_NN2(seq,c_seq=C_seq,dnac1=oligoConc, dnac2=oligoConc, Na=Na, Mg=Mg, Tris=tris,saltcorr=use_salt_correction, T_dG=T,
                       nn_table=NN_table, tmm_table=TMM_table, imm_table=IMM_table, de_table=DE_table,
                       dsCharge=dsCharge, F=f, disp_flag=debug, encerclingCoeff=encerclingCoeff
                       ,loop=loop, encercling=encercling_enable,rate=rate0,oligotype=oligotype , encerclingBpSalt=encerclingBpSalt
                             ,dss=ssDNA_d,bpss=ssDNA_bp,Sss=ssDNA_S,NNEnergyMod=modify,LAST_2bp=reduce_by_IN_on_last_2_bases,selfcomp=is_self_complementary(seq))
        local_salt_correction = use_salt_correction
        if use_salt_correction == 8:
            local_salt_correction = 5
        Bio_Tm = mt.Tm_NN(seq, dnac1=oligoConc, dnac2=oligoConc, Na=Na, Tris=tris, Mg=Mg, saltcorr=local_salt_correction,nn_table=NN_table, tmm_table=TMM_table, imm_table=IMM_table, de_table=DE_table,selfcomp=is_self_complementary(seq))
        if (Output == 'Toff'):
            if (debug) :
                print("Oligo-Sequence Toff (s)")
            if (debug == 0):
                print("%g\t%g Tm(BioSeq) %g" % (f,trep1[0][0,0],Bio_Tm));
            else :
                print("%s\t%g Tm(BioSeq) %g" % (seq,trep1[0][0,0],Bio_Tm));
        elif (Output == 'koff'):
            if (debug) :
                print("Oligo-Sequence koff (s-1)")
            print("%s\t%g" % (seq,trep1[1]));
        elif (Output == 'kon'):
            if (debug) :
                print("Oligo-Sequence kon (microM)")
            print("%s\t%g" % (seq,trep1[2]));
        elif (Output == 'dGF'):
            if (debug) :
                print("Oligo-Sequence dG(F) (kBT)")
            print("%s\t%g" % (seq,trep1[3]));
        elif (Output == 'Tm'):
            if (debug) :
                print("Oligo-Sequence Tm (C)\"\n")
            print("%s\t%g" % (seq,trep1[4]));

    elif (fmin >= 0) or (fmax >= 0) or (fstep >= 0):
        for f in numpy.arange( fmin, fmax, fstep):
            if (f != fmin):
                debug = 0
            trep1 = Tm_NN2(seq,c_seq=C_seq,dnac1=oligoConc, dnac2=oligoConc, Na=Na, Mg=Mg, Tris=tris, saltcorr=use_salt_correction, T_dG=T,
                           nn_table=NN_table, tmm_table=TMM_table, imm_table=IMM_table, de_table=DE_table,
                           dsCharge=dsCharge, F=f, disp_flag=debug, encerclingCoeff=encerclingCoeff , encerclingBpSalt=encerclingBpSalt
                           ,loop=loop, encercling=encercling_enable,rate=rate0,oligotype=oligotype,dss=ssDNA_d,bpss=ssDNA_bp,Sss=ssDNA_S,NNEnergyMod=modify,LAST_2bp=reduce_by_IN_on_last_2_bases,selfcomp=is_self_complementary(seq))
            if (Output == 'Toff'):            
                if (f == fmin):
                    print("\n-lx \"Force (pN)\"\n-ly \"Displacement time (s)\"\n")
                print("%g\t%g" % (f,trep1[0]));
            elif (Output == 'koff'):
                if (f == fmin):
                    print("\n-lx \"Force (pN)\"\n-ly \"Koff (s-1)\"\n")
                print("%g\t%g" % (f,trep1[1]));
            elif (Output == 'kon'):
                if (f == fmin):
                    print("\n-lx \"Force (pN)\"\n-ly \"Kon (microM)\"\n")
                print("%g\t%g" % (f,trep1[2]));
            elif (Output == 'dGF'):
                if (f == fmin):
                    print("\n-lx \"Force (pN)\"\n-ly \"dGF (kBT)\"\n")
                print("%g\t%g" % (f,trep1[3]));
            elif (Output == 'Tm'):
                if (f == fmin):
                    print("\n-lx \"Force (pN)\"\n-ly \"Tm (oC)\"\n")
                print("%g\t%g" % (f,trep1[4]));
        if kon == 1:
            debug = 0
            for f in numpy.arange( fmin, fmax, fstep):
                trep1 = Tm_NN2(seq,c_seq=C_seq,dnac1=oligoConc, dnac2=oligoConc, Na=Na, Mg=Mg, Tris=tris, saltcorr=use_salt_correction, T_dG=T,
                               nn_table=NN_table, tmm_table=TMM_table, imm_table=IMM_table, de_table=DE_table,
                               dsCharge=dsCharge, F=f, disp_flag=debug, encerclingCoeff=encerclingCoeff, encerclingBpSalt=encerclingBpSalt
                               ,loop=loop, encercling=encercling_enable,rate=rate0,oligotype=oligotype,dss=ssDNA_d,bpss=ssDNA_bp,Sss=ssDNA_S,NNEnergyMod=modify,LAST_2bp=reduce_by_IN_on_last_2_bases,selfcomp=is_self_complementary(seq))
                if (f == fmin):
                    print("\n-m 1\n-src \"Force (pN) koff (s-1)\"\n")
                print("%g\t%g" % (f,trep1[1]))
            for f in numpy.arange( fmin, fmax, fstep):
                trep1 = Tm_NN2(seq,c_seq=C_seq,dnac1=oligoConc, dnac2=oligoConc, Na=Na, Mg=Mg, Tris=tris, saltcorr=use_salt_correction, T_dG=T,
                               nn_table=NN_table, tmm_table=TMM_table, imm_table=IMM_table, de_table=DE_table,
                               dsCharge=dsCharge, F=f, disp_flag=debug, encerclingCoeff=encerclingCoeff, encerclingBpSalt=encerclingBpSalt
                               ,loop=loop, encercling=encercling_enable,rate=rate0,oligotype=oligotype,dss=ssDNA_d,bpss=ssDNA_bp,Sss=ssDNA_S,NNEnergyMod=modify,LAST_2bp=reduce_by_IN_on_last_2_bases,selfcomp=is_self_complementary(seq))
                if (f == fmin):
                    print("\n-m 1\n-src \"Force (pN) kon (x 10^6 M)\"\n")
                print("%g\t%g" % (f,trep1[2]))
            for f in numpy.arange( fmin, fmax, fstep):
                trep1 = Tm_NN2(seq,c_seq=C_seq,dnac1=oligoConc, dnac2=oligoConc, Na=Na, Mg=Mg, Tris=tris, saltcorr=use_salt_correction, T_dG=T,
                               nn_table=NN_table, tmm_table=TMM_table, imm_table=IMM_table, de_table=DE_table,
                               dsCharge=dsCharge, F=f, disp_flag=debug, encerclingCoeff=encerclingCoeff, encerclingBpSalt=encerclingBpSalt
                               ,loop=loop, encercling=encercling_enable,rate=rate0,oligotype=oligotype,dss=ssDNA_d,bpss=ssDNA_bp,Sss=ssDNA_S,NNEnergyMod=modify,LAST_2bp=reduce_by_IN_on_last_2_bases,selfcomp=is_self_complementary(seq))
                if (f == fmin):
                    print("\n-m 1\n-src \"Force (pN) DeltaG(F) (kBT)\"\n")
                print("%g\t%g" % (f,trep1[3]))


    elif (Tmin >= 0) or (Tmax >= 0) or (Tstep >= 0):
        for T in numpy.arange( Tmin, Tmax, Tstep):
            if (T != Tmin):
                debug = 0
            trep1 = Tm_NN2(seq,c_seq=C_seq,dnac1=oligoConc, dnac2=oligoConc, Na=Na, Mg=Mg, Tris=tris, saltcorr=use_salt_correction, T_dG=T,
                           nn_table=NN_table, tmm_table=TMM_table, imm_table=IMM_table, de_table=DE_table,
                           dsCharge=dsCharge, F=f, disp_flag=debug, encerclingCoeff=encerclingCoeff, encerclingBpSalt=encerclingBpSalt
                           ,loop=loop, encercling=encercling_enable,rate=rate0,oligotype=oligotype,dss=ssDNA_d,bpss=ssDNA_bp,Sss=ssDNA_S,NNEnergyMod=modify,LAST_2bp=reduce_by_IN_on_last_2_bases,selfcomp=is_self_complementary(seq))
            if (T == Tmin):
                print("\n-lx \"Temperature (C)\"\n-ly \"Displacement time (s)\"\n")
            print("%g\t%g" % (T,trep1[0]))

    return




if __name__ == "__main__":
    main()
