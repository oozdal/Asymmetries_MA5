from MyPySLHA import *
import pyslha
import fnmatch
import os, re

datafilepath     = os.path.abspath("/Users/oozdal/projects/SSM_WpZp/data/TEST.csv")
maindirpath = os.path.abspath("/Users/oozdal/projects/SSM_WpZp/")

numberofoutputdir = len(fnmatch.filter(os.listdir(maindirpath), 'run_*' ))

for i in range(1,numberofoutputdir+1): 
    outputdir = "run_"+"%.2d" % i
    print(outputdir)
    outputfilefullpath = outputdir +"/"+"run_"+"%.2d" % i +"_tag_1_banner.txt"
    print(outputfilefullpath)

    FullSysErrPath = outputdir + "/parton_systematics.log"

    SLHA = MyPySLHA()
    if SLHA.CheckLHAexist(outputfilefullpath) == True:
        SLHA.LoadLHAFile(outputfilefullpath)
        if    SLHA.Read_Decays(32,[-11,11]) == None:            SLHA.BRZptoee     = 0
        else : SLHA.BRZptoee     = SLHA.Read_Decays(32,[-11,11])
        if SLHA.Read_Decays(32,[-13,13]) == None:               SLHA.BRZptomumu   = 0
        else : SLHA.BRZptomumu   = SLHA.Read_Decays(32,[-13,13])
        SLHA.BRZptoll     = SLHA.BRZptoee + SLHA.BRZptomumu
        SLHA.MZp          = SLHA.ParticleMass
        SLHA.Zptotalwidth = SLHA.totalwidth

        SLHA.ReadxSection(outputfilefullpath)
        SLHA.ReadUncertainty(FullSysErrPath) 

        SLHA.MWp = SLHA.allcontent.blocks["MASS"][34]

        g= open(datafilepath,"a+")
        g.write(" %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E" % (SLHA.MZp, SLHA.xsection, SLHA.Zptotalwidth, SLHA.BRZptoee, SLHA.BRZptomumu, SLHA.BRZptoll, SLHA.upscale, SLHA.downscale, SLHA.upPDF, SLHA.downPDF, SLHA.MWp))
        g.write('\n')
        g.close()     
