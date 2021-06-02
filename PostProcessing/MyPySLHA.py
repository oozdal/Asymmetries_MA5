import pandas as pd
import numpy as np
import os
import shutil
import sys
# sys.path.append("/home/oo4g19/softwares/pyslha-3.2.3/")
import pyslha
from math import *
from scipy.interpolate import interp1d
import subprocess
#import xslha
import re
from decimal import Decimal


class MyPySLHA():
    def __init__(self):
        pass

    def LoadLesHouches(self, LesHouchesPath):
        self.LesHouches = xslha.read("LesHouches.in.sUMSSM_QH")

    def CheckLHAexist(self, LHAPath):
        return os.path.isfile(LHAPath)

    def LoadLHAFile(self, LHAPath):
        self.allcontent = pyslha.read(LHAPath, ignoreblocks=['SPINFO'])

    def WriteNewSLHAFile(self, newLHAPath, newLHAcontent, precision_val):
        self.newLHA = pyslha.writeSLHAFile(
            newLHAPath, newLHAcontent, precision=precision_val)
#        self.newLHA = pyslha.writeSLHA(newLHAcontent,ignorenobr=True)

    def Read_Block(self, BlockName, id1=None, id2=None, id3=None):
        self.BlockName = BlockName
        self.id1 = id1
        self.id2 = id2
        self.id3 = id3

        self.list_all_blocks = self.allcontent.blocks
        self.pyslha_Block = self.allcontent.blocks[self.BlockName]

        if self.id1 != None and self.id2 == None and self.id3 == None:
            self.VarValue1 = self.pyslha_Block[self.id1]
            return self.VarValue1

        elif self.id1 != None and self.id2 != None and self.id3 == None:
            self.VarValue2 = self.pyslha_Block[self.id1, self.id2]
            return self.VarValue2

        else:
            self.VarValue3 = self.pyslha_Block[self.id1, self.id2, self.id3]
            return self.VarValue3

    def Read_Decays(self, PySLHAparticleID, DecayProduct=None):

        self.PySLHAparticle = self.allcontent.decays[PySLHAparticleID]

        self.DecayProduct = DecayProduct

        self.list_decaysmodes = self.PySLHAparticle.decays
        self.totalwidth = self.PySLHAparticle.totalwidth
        self.ParticleMass = self.PySLHAparticle.mass

        for i in range(len(self.list_decaysmodes)):
            self.decaymode = self.list_decaysmodes[i]
            if self.decaymode.ids == self.DecayProduct:
                return self.decaymode.br

    def CheckMicrOMEGABlock(self):
        if self.allcontent.blocks.has_key("MICROMEGAS") == True:
            self.MicrOMEGASblock = self.allcontent.blocks["MICROMEGAS"]
            if self.MicrOMEGASblock.items() != []:
                return True
            else:
                return False
        else:
            return False

    def NormalizeDMResults(self):

        self.PlanckResult = 0.1187

        if self.CheckMicrOMEGABlock() == True:

            self.Relic_Density = self.Read_Block("MICROMEGAS", 700)
            self.RD_difference = abs(self.Relic_Density-self.PlanckResult)

            self.SI_proton = self.Read_Block("MICROMEGAS", 201)
            self.SI_neutron = self.Read_Block("MICROMEGAS", 203)

            self.SD_proton = self.Read_Block("MICROMEGAS", 202)
            self.SD_neutron = self.Read_Block("MICROMEGAS", 204)

            self.sigmaV = self.Read_Block("MICROMEGAS", 306)

            self.DirectDetec_NormFactor = min(
                1, self.Relic_Density/self.PlanckResult)
            self.IndirecDetec_NormFactor = min(
                1, (self.Relic_Density/self.PlanckResult)**2)

            self.NormSI_proton = self.SI_proton*self.DirectDetec_NormFactor
            self.NormSI_neutron = self.SI_neutron*self.DirectDetec_NormFactor

            self.SD_proton = self.SD_proton*self.DirectDetec_NormFactor
            self.SD_neutron = self.SD_neutron*self.DirectDetec_NormFactor

            self.Norm_sigmaV = self.sigmaV*self.IndirecDetec_NormFactor

    def MassesFromMicrOMEGAs(self):

        if self.CheckMicrOMEGABlock() == True:
            self.MICROMEGASblock = self.allcontent.blocks["MICROMEGAS"]

            self.MZp = self.Read_Block("MICROMEGAS", 401)
            self.MWp = self.Read_Block("MICROMEGAS", 402)

            self.Scot_el = self.Read_Block("MICROMEGAS", 403)
            self.Scot_mu = self.Read_Block("MICROMEGAS", 404)
            self.Scot_tau = self.Read_Block("MICROMEGAS", 405)

            self.MDD = self.Read_Block("MICROMEGAS", 406)
            self.MDS = self.Read_Block("MICROMEGAS", 407)
            self.MDB = self.Read_Block("MICROMEGAS", 408)

            self.mh1 = self.Read_Block("MICROMEGAS", 409)
            self.mh2 = self.Read_Block("MICROMEGAS", 410)
            self.mh3 = self.Read_Block("MICROMEGAS", 411)

            self.mA1 = self.Read_Block("MICROMEGAS", 412)
            self.mA2 = self.Read_Block("MICROMEGAS", 413)

            self.mhp = self.Read_Block("MICROMEGAS", 414)
            self.mhm = self.Read_Block("MICROMEGAS", 415)

            if self.MICROMEGASblock.has_key(900) == True:
                self.pptoWpWp = self.Read_Block("MICROMEGAS", 900)
            else:
                self.pptoWpWp = None

            if self.MICROMEGASblock.has_key(901) == True:
                self.pptoWpdd = self.Read_Block("MICROMEGAS", 901)
            else:
                self.pptoWpdd = None

            if self.MICROMEGASblock.has_key(902) == True:
                self.pptodddd = self.Read_Block("MICROMEGAS", 902)
            else:
                self.pptodddd = None

            if self.MICROMEGASblock.has_key(903) == True:
                self.pptoZp = self.Read_Block("MICROMEGAS", 903)
            else:
                self.pptoZp = None

    def MassConstraints(self):

        self.MZp_Bound = False
        self.Mh2_Bound = False
        self.Mh3_Bound = False
        self.fshm2_Bound = False
        self.MDD_Bound = False
        self.Scot_DM_Bound = False
        self.Scot_Hierarchy = False
        self.dprime_Hierarchy = False
        self.DMmass_Bound = False
        self.WpMass_Bound = False

        if self.CheckMicrOMEGABlock() == True:
            if self.Read_Block("MICROMEGAS", 401) > 3750.0 and self.Read_Block("MICROMEGAS", 401) < 5500.0:
                self.MZp_Bound = True
            if self.Read_Block("MICROMEGAS", 403) < 2000:
                self.Scot_DM_Bound = True
            if (self.Read_Block("MICROMEGAS", 403) < self.Read_Block("MICROMEGAS", 404) < self.Read_Block("MICROMEGAS", 405)):
                self.Scot_Hierarchy = True
            if (self.Read_Block("MICROMEGAS", 406) < self.Read_Block("MICROMEGAS", 407) < self.Read_Block("MICROMEGAS", 408)):
                self.dprime_Hierarchy = True
#            if self.Read_Block("MICROMEGAS",410) < 500.0:
#                self.Mh2_Bound = True
#            if self.Read_Block("MICROMEGAS",411) < 250000.:
#                self.Mh3_Bound = True
#            if self.Read_Block("MICROMEGAS",415) > self.Read_Block("MASS",25)/2.:
#                self.fshm2_Bound = True
#            if self.Read_Block("MICROMEGAS",406) < 250000.:
#                self.MDD_Bound = True

            self.AssignDMMass()
            if self.DMmass < 2000:
                self.DMmass_Bound = True
            if self.MWp > self.DMmass:
                self.WpMass_Bound = True

#        self.CheckConstraints = self.MZp_Bound and self.fshm2_Bound and self.Mh3_Bound and self.MDD_Bound
        self.CheckConstraints = self.MZp_Bound and self.DMmass_Bound and self.WpMass_Bound
        return self.CheckConstraints

    def FindDMParticle(self):
        if (self.Read_Block("MICROMEGAS", 403) < self.Read_Block("MICROMEGAS", 404)) and (self.Read_Block("MICROMEGAS", 403) < self.Read_Block("MICROMEGAS", 405)):
            self.DMParticle = "~ne"
            self.AntiDMParticle = "~ne~"
        elif (self.Read_Block("MICROMEGAS", 404) < self.Read_Block("MICROMEGAS", 403)) and (self.Read_Block("MICROMEGAS", 404) < self.Read_Block("MICROMEGAS", 405)):
            self.DMParticle = "~nm"
            self.AntiDMParticle = "~nm~"
        elif (self.Read_Block("MICROMEGAS", 405) < self.Read_Block("MICROMEGAS", 403)) and (self.Read_Block("MICROMEGAS", 405) < self.Read_Block("MICROMEGAS", 404)):
            self.DMParticle = "~nt"
            self.AntiDMParticle = "~nt~"
        else:
            self.DMParticle = "Unknown"
            self.AntiDMParticle = "Unknown"

    def AssignDMMass(self):
        self.FindDMParticle()
        if self.DMParticle == "~ne":
            self.DMmass = self.Read_Block("MICROMEGAS", 403)
        elif self.DMParticle == "~nm":
            self.DMmass = self.Read_Block("MICROMEGAS", 404)
        elif self.DMParticle == "~nt":
            self.DMmass = self.Read_Block("MICROMEGAS", 405)
        else:
            self.DMmass = None

    def RelicDensity_Constraint(self):

        if self.CheckMicrOMEGABlock() == True:
            self.Relic_Density = self.Read_Block("MICROMEGAS", 700)
            if self.Relic_Density <= 0.2:
                self.Relic_Density_Bound = True
            else:
                self.Relic_Density_Bound = False

        else:
            self.Relic_Density_Bound = False

        return self.Relic_Density_Bound

    def CalErrorPercent(self, experimental_val, theoretical_val):
        self.experimental_val = experimental_val
        self.theoretical_val = theoretical_val

        self.ErrorPercent = abs(
            (self.experimental_val-self.theoretical_val)/self.theoretical_val)*100

        return self.ErrorPercent

    def ReadMicrOmegasOutput(self, MicrOmegasFileName):
        self.MicrOmegasFileName = MicrOmegasFileName
        self.LoadMicrOmegasFile = pd.read_csv(self.MicrOmegasFileName, skipinitialspace=True, header=1, sep="   | # ", names=[
                                              "Coding", "Values", "Comment"], engine='python')

        for index in range(len(self.LoadMicrOmegasFile.index)):
            if self.LoadMicrOmegasFile["Coding"][index] == 700:
                self.RelicDensity = "{:.6E}".format(
                    float(self.LoadMicrOmegasFile["Values"][0]))
                self.RelicDensity = float(self.RelicDensity)
            if self.LoadMicrOmegasFile["Coding"][index] == 903:
                self.xsecpptoZp = "{:.6E}".format(
                    float(self.LoadMicrOmegasFile["Values"][index]))
                self.xsecpptoZp = float(self.xsecpptoZp)
            if self.LoadMicrOmegasFile["Coding"][index] == 201:
                self.SIproton = "{:.6E}".format(
                    float(self.LoadMicrOmegasFile["Values"][index]))
                self.SIproton = float(self.SIproton)
            if self.LoadMicrOmegasFile["Coding"][index] == 203:
                self.SIneutron = "{:.6E}".format(
                    float(self.LoadMicrOmegasFile["Values"][index]))
                self.SIneutron = float(self.SIneutron)
            if self.LoadMicrOmegasFile["Coding"][index] == 306:
                self.sigmaV = "{:.6E}".format(
                    float(self.LoadMicrOmegasFile["Values"][index]))
                self.sigmaV = float(self.sigmaV)
            if self.LoadMicrOmegasFile["Coding"][index] == 904:
                self.xsecCha1Cha1 = "{:.6E}".format(
                    float(self.LoadMicrOmegasFile["Values"][index]))
                self.xsecCha1Cha1 = float(self.xsecCha1Cha1)
            if self.LoadMicrOmegasFile["Coding"][index] == 905:
                self.xsecChi2Cha1 = "{:.6E}".format(
                    float(self.LoadMicrOmegasFile["Values"][index]))
                self.xsecChi2Cha1 = float(self.xsecChi2Cha1)
            if self.LoadMicrOmegasFile["Coding"][index] == 906:
                self.xsecChi2AntiCha1 = "{:.6E}".format(
                    float(self.LoadMicrOmegasFile["Values"][index]))
                self.xsecChi2AntiCha1 = float(self.xsecChi2Cha1)

    def ReadChannels(self, ChannelFileName):
        self.ChannelFileName = ChannelFileName
        self.LoadChannelFile = pd.read_csv(
            self.ChannelFileName, skipinitialspace=True, header=1, sep="%", names=["Percent", "Channels"])
        self.LoadChannelFile = self.LoadChannelFile.set_index(
            ["Percent", "Channels"])

    def CheckHiggsFunnel(self):
        HiggsFunnelExist = "~ne ~ne~ ->H0 H0 " in self.LoadChannelFile.index.levels[1]
        if HiggsFunnelExist == True:
            self.ChannelswithIndex = self.LoadChannelFile.index
            for index in range(len(self.ChannelswithIndex)):
                if self.ChannelswithIndex[index][1] == "~ne ~ne~ ->H0 H0 ":
                    self.HiggsFunnelPercentage = self.ChannelswithIndex[index][0]
            return self.HiggsFunnelPercentage
        else:
            self.HiggsFunnelPercentage = None
            return self.HiggsFunnelPercentage

    def CheckChannelContribution(self, ChannelName):
        self.ChannelName = ChannelName

        ChannelExist = self.ChannelName in self.LoadChannelFile.index.levels[1]
        if ChannelExist == True:
            self.ChannelswithIndex = self.LoadChannelFile.index
            for index in range(len(self.ChannelswithIndex)):
                if self.ChannelswithIndex[index][1] == self.ChannelName:
                    self.ChannelPercentage = self.ChannelswithIndex[index][0]
            return self.ChannelPercentage
        else:
            self.ChannelPercentage = 0
            return self.ChannelPercentage

    def Erase(self, FilePathToErase):
        self.FilePathToErase = FilePathToErase
        if os.path.isfile(self.FilePathToErase) == True:
            os.remove(self.FilePathToErase)

    def RD_difference(self):
        self.PlanckResult = 0.1187
        self.diff_RD = abs(self.Relic_Density-self.PlanckResult)

        return self.diff_RD

    def Check_RD_diff(self, RD_diff_list):
        self.RD_diff_list = RD_diff_list
        if len(self.RD_diff_list) >= 2:
            if self.RD_diff_list[-1] <= self.RD_diff_list[-2]:
                self.Logic_RD_diff = True
                self.RD_diff_list.pop(0)
            else:
                self.Logic_RD_diff = False
        elif len(self.RD_diff_list) == 1:
            self.Logic_RD_diff = True

#        return self.Logic_RD_diff
        return True

    def CheckXENON1TLimit(self):

        self.LoadExpData(
            "XENON1T", "/home/phylab/hepwork/ExpDATA/expdata/XENON1T.dat")

        self.nearest_ObservedLim_xsection = []
        self.nearest_ObservedLim_DMmass = []

        self.nearest_ObservedLim_xsection.append(self.XENON1Tdata["XENON1TXSECTION"][self.Closest(
            self.XENON1Tdata["XENON1TWIMPMASS"], self.DMmass)])
        self.nearest_ObservedLim_DMmass.append(self.XENON1Tdata["XENON1TWIMPMASS"][self.Closest(
            self.XENON1Tdata["XENON1TWIMPMASS"], self.DMmass)])

        if ((self.nearest_ObservedLim_xsection[0] >= self.SI_proton) and (self.nearest_ObservedLim_xsection[0] >= self.SI_neutron)):
            return True
        else:
            return False

    def ZpMassLimit(self):
        #################################################
        #        self.given_gR = self.gR

        #        self.gR_values_given = [0.37, self.gL]
        #        self.Zp_values_given = [5000., 3950.]

        #        self.y1 = interp1d(self.gR_values_given, self.Zp_values_given)
        #        self.Zp_limit = self.y1(self.given_gR)
        #        self.Zp_limit = self.Zp_limit.tolist()
        ##################################################

        self.LoadExpData("observedxsection13TeVZpll",
                         "/home/phylab/hepwork/ExpDATA/expdata/observedxsection13TeVZpll.dat")
        self.LoadExpData("expectedxsection13TeVZpll",
                         "/home/phylab/hepwork/ExpDATA/expdata/expectedxsection13TeVZpll.dat")

        self.nearest_ObservedLim_xsection = []
        self.nearest_ExpectedLim_xsection = []

        self.nearest_ObservedLim_Zpmass = []
        self.nearest_ExpectedLim_Zpmass = []

        self.MZp_TeV = self.MZp/1000.

        self.nearest_ObservedLim_xsection.append(self.Zptoll_ObservedLimData["xsection"][self.Closest(
            self.Zptoll_ObservedLimData["ZRmass"], self.MZp_TeV)])
        self.nearest_ExpectedLim_xsection.append(self.Zptoll_ExpectedLimData["xsection"][self.Closest(
            self.Zptoll_ExpectedLimData["ZRmass"], self.MZp_TeV)])

        self.nearest_ObservedLim_Zpmass.append(self.Zptoll_ObservedLimData["ZRmass"][self.Closest(
            self.Zptoll_ObservedLimData["ZRmass"], self.MZp_TeV)])
        self.nearest_ExpectedLim_Zpmass.append(self.Zptoll_ExpectedLimData["ZRmass"][self.Closest(
            self.Zptoll_ExpectedLimData["ZRmass"], self.MZp_TeV)])

    def CheckZpMassLimit(self):
        ###########################################################
        #        self.ZpMassLimit()

        #        if self.CheckMicrOMEGABlock() == True:
        #            if self.MZp < self.ZpMassLimit():
        #                self.LogicZpMassLimit = False
        #            elif self.MZp >= self.ZpMassLimit():
        #                self.LogicZpMassLimit = True
        #        else:
        #            self.LogicZpMassLimit = False

        #        return self.LogicZpMassLimit
        ###########################################################

        self.ZpMassLimit()

        self.BRZptoll = self.ReadDecayChannel(
            " Zp -> e-,e+") + self.ReadDecayChannel(" Zp -> mu-,mu+")

        if self.nearest_ObservedLim_xsection[0] >= self.pptoZp*self.BRZptoll:
            return True
        else:
            return False

#        return True

    def RunHiggsBounds(self, InputFile):
        h = subprocess.call(
            "timeout 20 ./../../../hepwork/HiggsBounds-4.3.1/HiggsBounds LandH SLHA 3 1 "+str(InputFile), shell=True)

    def CheckHiggsBounds(self):
        if self.Read_Block("HIGGSBOUNDSRESULTS", 1, 2) == 1:
            self.HiggsBounds_Logic = True
        else:
            self.HiggsBounds_Logic = False

    def RunHiggsSignals(self, InputFile):
        g = subprocess.call(
            "timeout 20 ./../../../hepwork/HiggsSignals-1.4.0/HiggsSignals latestresults mass 2 SLHA 3 1 "+str(InputFile), shell=True)

    def ReadDecayFile(self, DecayFileName, FlagNo):
        self.DecayFileName = DecayFileName
        self.FlagNo = FlagNo

        if self.FlagNo == 1:
            self.DecayFile = pd.read_csv(self.DecayFileName, skipinitialspace=False,
                                         sep="->|#", names=["BRval", "BRChannels"], engine='python')

            for i in range(len(self.DecayFile)):
                self.temp_str1 = self.DecayFile["BRval"][i]
                self.temp_str2 = self.DecayFile["BRChannels"][i]

                self.temp_str1_split = self.temp_str1.split()

                if len(self.temp_str1_split) >= 2:
                    self.DecayFile["BRval"][i] = self.temp_str1_split[0]
                    self.DecayFile["BRChannels"][i] = self.temp_str1_split[1] + \
                        " ->"+self.temp_str2

            self.BranchingRatios = self.DecayFile["BRChannels"]

        if self.FlagNo == 2:
            self.DecayFile = pd.read_csv(self.DecayFileName, sep="  |   ", comment="Branching", names=[
                                         "BRvalues", "PartialWidths", "ChannelNames"], engine='python')
            self.ChannelNames = self.DecayFile["ChannelNames"]
            self.PartialWidths = self.DecayFile["PartialWidths"]
            self.BRvalues = self.DecayFile["BRvalues"]

    def ReadDecayChannel(self, CheckChannelName=None):

        self.CheckChannelName = CheckChannelName
        if self.CheckChannelName in self.ChannelNames.unique():
            for i in range(len(self.DecayFile)):
                if self.DecayFile["ChannelNames"][i] == self.CheckChannelName:
                    self.BRval = float(self.DecayFile["BRvalues"][i])
                    self.PartialWidthWithUnit = self.DecayFile["PartialWidths"][i]
                    self.PartialWidth_temp = self.PartialWidthWithUnit.split()
                    self.PartialWidth = self.PartialWidth_temp[0]
                    return self.BRval
        else:
            self.BRval = 0
            return self.BRval

    def ReadBRthroughDecayFile(self, BRchannel=None):
        self.BRchannel = BRchannel
        if self.BRchannel in self.BranchingRatios.unique():
            for i in range(len(self.DecayFile)):
                if self.DecayFile["BRChannels"][i] == self.BRchannel:
                    self.BRval = float(self.DecayFile["BRval"][i])
                    return self.BRval
        else:
            self.BRval = 0
            return self.BRval

    def DecayRatesThroughDecayFile(self):
        self.BRHiggstogammagamma = self.ReadBRthroughDecayFile("H0 -> A,A")
        self.DRwithBlocks_Higgstogammagamma = self.mh0_totalwidth*self.BRHiggstogammagamma

        self.BRHiggstoGluGlu = self.ReadBRthroughDecayFile("H0 -> G,G")
        self.DRwithBlocks_HiggstoGluGlu = self.mh0_totalwidth*self.BRHiggstoGluGlu

        self.DR_HiggstoGluGlu_SMvalue = 3.5e-4
#        self.Kgluglu_Rate = self.DRwithBlocks_HiggstoGluGlu/self.DR_HiggstoGluGlu_SMvalue
        self.Kgluglu_Rate = self.DRwithBlocks_HiggstoGluGlu/self.DR_HiggstoGluGlu_SM
        self.Kgammagamma_Rate = self.DRwithBlocks_Higgstogammagamma / \
            self.DR_HiggstoGammaGamma_SM
        self.Mgammagamma_Rate = self.Kgluglu_Rate*self.Kgammagamma_Rate

    def Convertcm2topb(self, xsectioncm2):
        xsectionpb = xsectioncm2*1e36
        return xsectionpb

    def LoadExpData(self, ExperimentName, ExpDatapath):
        if ExperimentName == "XENON1T":
            self.XENON1Tdata = pd.read_csv(ExpDatapath, sep=",", header=None, names=[
                                           "XENON1TWIMPMASS", "XENON1TXSECTION"])
            self.XENON1Tdata["XENON1TXSECTION"] = self.Convertcm2topb(
                self.XENON1Tdata["XENON1TXSECTION"])
        elif ExperimentName == "CMS_WRtoqqee_ObservedLimit":
            self.WRtoqqee_ObservedLimData = pd.read_csv(ExpDatapath, sep=",", header=None, names=[
                                                        "WRmass", "xsection"])  # This data is in fb unit
        elif ExperimentName == "CMS_WRtoqqee_ExpectedLimit":
            self.WRtoqqee_ExpectedLimData = pd.read_csv(ExpDatapath, sep=",", header=None, names=[
                                                        "WRmass", "xsection"])  # This data is in fb unit
        elif ExperimentName == "CMS_WRtoqqmumu_ObservedLimit":
            self.WRtoqqmumu_ObservedLimData = pd.read_csv(ExpDatapath, sep=",", header=None, names=[
                                                          "WRmass", "xsection"])  # This data is in fb unit
        elif ExperimentName == "CMS_WRtoqqmumu_ExpectedLimit":
            self.WRtoqqmumu_ExpectedLimData = pd.read_csv(ExpDatapath, sep=",", header=None, names=[
                                                          "WRmass", "xsection"])  # This data is in fb unit
        elif ExperimentName == "ATLAS_Efficiency_SS_ee.csv":
            self.Efficiency_SS_ee = pd.read_csv(ExpDatapath, sep=",", header=11, names=[
                                                "WRmass", "NRmass", "Efficiency"])  # Masses are in TeV unit
        elif ExperimentName == "ATLAS_Efficiency_OS_ee.csv":
            self.Efficiency_OS_ee = pd.read_csv(ExpDatapath, sep=",", header=11, names=[
                                                "WRmass", "NRmass", "Efficiency"])  # Masses are in TeV unit
        elif ExperimentName == "ATLAS_Efficiency_SS_mumu.csv":
            self.Efficiency_SS_mumu = pd.read_csv(ExpDatapath, sep=",", header=11, names=[
                                                  "WRmass", "NRmass", "Efficiency"])  # Masses are in TeV unit
        elif ExperimentName == "ATLAS_Efficiency_OS_mumu.csv":
            self.Efficiency_OS_mumu = pd.read_csv(ExpDatapath, sep=",", header=11, names=[
                                                  "WRmass", "NRmass", "Efficiency"])  # Masses are in TeV unit
        elif ExperimentName == "observedxsection13TeVZpll":
            self.Zptoll_ObservedLimData = pd.read_csv(ExpDatapath, sep=",", header=11, names=[
                                                      "ZRmass", "xsection"])  # Masses are in TeV unit
        elif ExperimentName == "expectedxsection13TeVZpll":
            self.Zptoll_ExpectedLimData = pd.read_csv(ExpDatapath, sep=",", header=11, names=[
                                                      "ZRmass", "xsection"])  # Masses are in TeV unit
        else:
            pass

    def Closest(self, list, Number):
        aux = []
        for valor in list:
            aux.append(abs(Number-valor))
        return aux.index(min(aux))

    def U1ChargeGenerator(self, alpha=None):
        self.alpha = alpha

	# Basis: Qu1, Qu2, Qu3, Qd1, Qd2, Qd3, Qe1, Qe2, Qe3, Qq1, Qq2, Qq3, Ql1, Ql2, Ql3, QHu, QHd, Qs
        self.ChargeSets = pd.read_csv('/Users/oozdal/projects/MagMoment/Scan/NonUMSSMChargeSets.csv', header=0,sep=", ", engine='python')

	# Excluding Qu2, Qu3, Qd2, Qd3, Qq2 and Qq3 since Qq1 = Qq2 = Qq3, Qu1 = Qu2 = Qu3 and Qd1 = Qd2 = Qd3
        self.ChargeSets = self.ChargeSets[['Qu1', 'Qd1','Qe1','Qe2','Qe3','Qq1','Ql1','Ql2','Ql3','QHu','QHd','Qs']] 
        self.ChargeSets.columns =['Qu', 'Qd','Qe1','Qe2','Qe3','Qq','Ql1','Ql2','Ql3','QHu','QHd','Qs']  # Renaming the column names

        n = np.array([i for i in range(56)])
	randcho = np.random.choice(n)

        if self.alpha == None: 

        	self.QQ_charge  = self.ChargeSets['Qq'][randcho] 
        	self.Ql1_charge = self.ChargeSets['Ql1'][randcho]   
        	self.Ql2_charge = self.ChargeSets['Ql2'][randcho]    
        	self.Ql3_charge = self.ChargeSets['Ql3'][randcho]  
        	self.QHu_charge = self.ChargeSets['QHu'][randcho]    
        	self.QHd_charge = self.ChargeSets['QHd'][randcho]    
        	self.Qd_charge  = self.ChargeSets['Qd'][randcho]  
        	self.Qu_charge  = self.ChargeSets['Qu'][randcho]
        	self.Qe1_charge = self.ChargeSets['Qe1'][randcho]   
        	self.Qe2_charge = self.ChargeSets['Qe2'][randcho]               
        	self.Qe3_charge = self.ChargeSets['Qe3'][randcho]                   
        	self.Qs_charge  = self.ChargeSets['Qs'][randcho]  

        if self.alpha != None:
		self.QQ_charge  = self.alpha/9.
		self.Ql1_charge = -3*self.alpha
		self.Ql2_charge = self.alpha
		self.Ql3_charge = -1*self.alpha
		self.QHu_charge = self.alpha
		self.QHd_charge = self.alpha
		self.Qd_charge  = (8*self.alpha)/9.
		self.Qu_charge  = (-10*self.alpha)/9.
		self.Qe1_charge = 4*self.alpha
		self.Qe2_charge = 0.
		self.Qe3_charge = 0.
		self.Qs_charge  = -2*self.alpha
    def ReadxSection(self, xsectionPath):
        file = open(xsectionPath, "r")
        for line in file:
            if line[3:20] == "Integrated weight":
                self.xsection = float(line[35:52])

    def ReadUncertainty(self, SysErrorPath):
        sysfile = open(SysErrorPath, "r")
        for line in sysfile:
            if line[6:21] == "scale variation":
                scalevar = re.split("\+|%|-|%", line)
                self.upscale = float(scalevar[1])
                self.downscale = float(scalevar[3])
            
            if line[2:15] == "PDF variation":
                pdfvar = re.split("\+|%|-|%", line)
                self.upPDF   = float(pdfvar[1]) 
                self.downPDF = float(pdfvar[3])
          

    def RenameAndCopy(self, currentLHAFilePath, newLHAFilePath):
        self.currentLHAFullPath = os.path.abspath(currentLHAFilePath)
        self.currentLHAdirPath = os.path.dirname(self.currentLHAFullPath)
        self.currentLHAFileName = os.path.basename(self.currentLHAFullPath)

        self.newLHAFullPath = os.path.abspath(newLHAFilePath)
        self.newLHAdirPath = os.path.dirname(self.newLHAFullPath)
        self.newLHAFileName = os.path.basename(self.newLHAFullPath)

        self.old_file = os.path.join(
            self.currentLHAdirPath, self.currentLHAFileName)
        self.new_file = os.path.join(self.newLHAdirPath, self.newLHAFileName)
        shutil.copy2(self.old_file, self.new_file)

    def RunSPheno(self, LHAInput):

        #        p = subprocess.call("./../../../softwares/SPheno-4.0.4/bin/SPhenosec_UMSSM "+str(LHAInput),shell=True)
        #        p = subprocess.call("./../../../../../home/oo4g19/softwares/SPheno-4.0.4/bin/SPhenosec_UMSSM "+str(LHAInput),shell=True)
        p = subprocess.call(
            "./../../../packages/SPheno-4.0.4/bin/SPhenoNonUMSSM "+str(LHAInput), shell=True)

    def MassBounds(self):

        self.SMlikeHiggsMass = False
        self.GluinoBound = False
        self.Chi1Bound = False
        self.Cha1Bound = False
        self.ZpMassBound = False
        self.ZZmixingBound = False

        self.mh1 = self.allcontent.blocks["MASS"][25]
        self.mh2 = self.allcontent.blocks["MASS"][35]
        self.mh3 = self.allcontent.blocks["MASS"][45]
        self.MZp = self.allcontent.blocks["MASS"][31]
        self.Gluino = self.allcontent.blocks["MASS"][1000021]
        self.Chi1 = abs(self.allcontent.blocks["MASS"][1000022])
        self.Cha1 = self.allcontent.blocks["MASS"][1000024]
        self.ZZmixing = abs(self.allcontent.blocks["ANGLES"][10])

        if (self.mh1 >= 122 and self.mh1 <= 128) or (self.mh2 >= 122 and self.mh2 <= 128) or (self.mh3 >= 122 and self.mh3 <= 128):
            self.SMlikeHiggsMass = True
        if (self.Gluino >= 1800):
            self.GluinoBound = True

        if (self.Cha1 >= 103.5):
            self.Cha1Bound = True

        if (self.MZp > 1000) and (self.MZp < 6000)	:
            self.ZpMassBound = True
        if (self.ZZmixing <= 1e-3):
            self.ZZmixingBound = True

        # and self.ZZmixingBound
        self.CheckConstraints = self.SMlikeHiggsMass and self.GluinoBound and self.Cha1Bound and self.ZpMassBound
        return self.CheckConstraints

    def G2Bound(self):
        self.G2Sign = False
        self.DAMU3sigma = False
        self.DAMU2sigma = False
        self.DAMU1sigma = False

        self.DAEL         = self.allcontent.blocks["SPHENOLOWENERGY"][20]
        self.DAMU         = self.allcontent.blocks["SPHENOLOWENERGY"][21]

        if self.DAEL < 0 and self.DAMU > 0: self.G2Sign = True

        if self.DAMU*10**10 > 3.4  and self.DAMU*10**10 < 55.6: self.DAMU3sigma = True
        if self.DAMU*10**10 > 12.7 and self.DAMU*10**10 < 44.7: self.DAMU2sigma = True
        if self.DAMU*10**10 > 20.7 and self.DAMU*10**10 < 36.7: self.DAMU1sigma = True

        if (self.G2Sign and self.DAMU3sigma):  return True
        else: return False        

    def CheckWhoIsLSP(self):
        self.LSPcontraint = False

#        if (self.allcontent.blocks["LSP"][1] == 1000022) or (self.allcontent.blocks["LSP"][1] == 1000012): self.LSPcontraint = True
        if (self.allcontent.blocks["LSP"][1] == 1000022):
            self.LSPcontraint = True
        else:
            self.LSPcontraint == False

        return self.LSPcontraint

    def CheckBPhysics(self):
        self.BtoXsgamma = False
        self.Bstomumu = False
        self.RBtotaunu = False
        self.BPhysics = False

        if (self.allcontent.blocks["FLAVORKITQFV"][200] >= 2.99E-4) and (self.allcontent.blocks["FLAVORKITQFV"][200] <= 3.87E-4):
            self.BtoXsgamma = True
        if (self.allcontent.blocks["FLAVORKITQFV"][4006] >= 1.1E-9) and (self.allcontent.blocks["FLAVORKITQFV"][4006] <= 6.4E-9):
            self.Bstomumu = True
        if (self.allcontent.blocks["FLAVORKITQFV"][503] >= 0.15) and (self.allcontent.blocks["FLAVORKITQFV"][503] <= 2.41):
            self.RBtotaunu = True
        self.BPhysics = self.BtoXsgamma and self.Bstomumu and self.RBtotaunu
        return self.BPhysics

    def LSPmass(self):

        self.LSPmass = []
        for i in range(len(self.data)):
            if self.param["LSP"][i] == 1000022:
                self.LSPmass.append(self.param["mchi1"][i])
            elif self.param["LSP"][i] == 1000012:
                self.LSPmass.append(self.param["MassSv1"][i])

        self.LSPmass = abs(np.array(self.LSPmass))

    def NLSPmass(self):

        self.NLSPmass = []
        for i in range(len(self.data)):
            if self.param["NLSP"][i] == 1000012:
                self.NLSPmass.append(self.param["MassSv1"][i])
            elif self.param["NLSP"][i] == 1000023:
                self.NLSPmass.append(self.param["mchi2"][i])
            elif self.param["NLSP"][i] == 1000024:
                self.NLSPmass.append(self.param["mcha1"][i])
            elif self.param["NLSP"][i] == 1000011:
                self.NLSPmass.append(self.param["stau1"][i])
            elif self.param["NLSP"][i] == 1000021:
                self.NLSPmass.append(self.param["gluino"][i])
            elif self.param["NLSP"][i] == 1000022:
                self.NLSPmass.append(self.param["mchi1"][i])
            else:
                self.NLSPmass.append(1e40)

    def RunMicrOMEGAs(self, LHAInput):

        #        p = subprocess.call("./../../../softwares/micromegas_5.0.8/secluded_UMSSM/myOmega "+str(LHAInput),shell=True)
        p = subprocess.call("./myOmega "+str(LHAInput), shell=True)

    def LHAwithDM(self, PathLHAfile, PathMicrOMEGAsResult, PathLHAfileWithDM):
        #        c = subprocess.call("timeout 5 cat "+str(PathLHAfile)+" "+str(PathMicrOMEGAsResult)+" > "+str(PathLHAfileWithDM),shell=True)

        filenames = [PathLHAfile, PathMicrOMEGAsResult]
        with open(PathLHAfileWithDM, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def chi2STU(self, sigma=2.):  # 1508.01671 p7 eq 10-12
        S = self.allcontent.blocks['SPHENOLOWENERGY'][2]
        T = self.allcontent.blocks['SPHENOLOWENERGY'][1]
        U = self.allcontent.blocks['SPHENOLOWENERGY'][3]
        from numpy import array, dot, sqrt
        from numpy.linalg import inv
        from scipy import stats, special
        x = array([[S-.05], [T-.09], [U-.01]])
        xT = x.transpose()
        Cij = array([[.0121, .0129, -.0071],
                     [.0129, .0169, -.0119],
                     [-.0071, -.0119, .0121]])
        invC = inv(Cij)
        # chi2 = xT.Cij-1.x
        xTinvCij = dot(xT, invC)
        self.chi2 = dot(xTinvCij, x)
        if float(self.chi2) <= stats.chi2.isf((1-special.erf(sigma/sqrt(2.))), 3):
            return True
        else:
            return False

    def Cha1Chi1Excl(self, test_Cha1, test_Chi1):
        self.test_Cha1 = abs(test_Cha1)
        self.test_Chi1 = abs(test_Chi1)
        exec(open("/home/oo4g19/softwares/ExpDATA/ExDATA.py").read())
        self.coord = []
        for i in range(len(observedATLAS_CharNeut_139fb_Cha1)):
            self.coord.append(
                (observedATLAS_CharNeut_139fb_Cha1[i], observedATLAS_CharNeut_139fb_Chi1[i]))
        poly = Polygon(self.coord)
        self.test_point = Point(self.test_Cha1, self.test_Chi1)
        if self.test_point.within(poly) == True:
            return True
        if self.test_point.within(poly) == False:
            return False

    def Chi1Content(self):
        self.NMIX11 = self.allcontent.blocks["NMNMIX"][1, 1]
        self.NMIX12 = self.allcontent.blocks["NMNMIX"][1, 2]
        self.NMIX13 = self.allcontent.blocks["NMNMIX"][1, 3]
        self.NMIX14 = self.allcontent.blocks["NMNMIX"][1, 4]
        self.NMIX15 = self.allcontent.blocks["NMNMIX"][1, 5]
        self.NMIX16 = self.allcontent.blocks["NMNMIX"][1, 6]
#        self.NMIX17 = self.allcontent.blocks["NMNMIX"][1, 7]
#        self.NMIX18 = self.allcontent.blocks["NMNMIX"][1, 8]
#        self.NMIX19 = self.allcontent.blocks["NMNMIX"][1, 9]

    def AddXSecBlock(self, SLHAFile, Destination, Particle1, Particle2, XSEC):
        self.SLHAFile = os.path.abspath(SLHAFile)
        self.Destination = Destination
        self.Particle1 = Particle1
        self.Particle2 = Particle2
        self.XSEC = XSEC

        self.xsecLine1 = "XSECTION  1.30E+04  2212 2212 2 "+self.Particle1 + \
            " "+self.Particle2+" # 10000 events, [pb], pythia8 for LO"
        self.xsecline2 = "  0  0  0  0  0  0    "+"%.8E" % self.XSEC+" SModelSv1.2.2"

        self.RenameAndCopy(self.SLHAFile, self.Destination)

        # Open a file with access mode 'a'
        file_object = open(self.Destination, 'a')

        # Append 'hello' at the end of file
        file_object.write("\n")
        file_object.write(self.xsecLine1+"\n")
        file_object.write(self.xsecline2)
        file_object.write("\n")

        # Close the file
        file_object.close()

    def AssignValues(self, LesHouchesFileFullPath, UpdatedLesHouchesFileFullPath):
        self.LesHouchesFileFullPath = LesHouchesFileFullPath
        self.UpdatedLesHouchesFileFullPath = UpdatedLesHouchesFileFullPath

        if self.CheckFile(self.LesHouchesFileFullPath) == True:
            with open(self.LesHouchesFileFullPath, 'r') as self.file:
                self.filedata = self.file.read()

            # Replace the target string
            for i in range(self.numberofparameters):
                self.filedata = self.filedata.replace(self.freeparamsdata["ParamName"][i], str(
                    "{:.6E}".format(self.params[str(self.freeparamsdata["ParamName"][i])])))

            # Write the file out again
            with open(self.UpdatedLesHouchesFileFullPath, 'w') as updatedfile:
                self.filedata = updatedfile.write(self.filedata)
            updatedfile.close()
        else:
            print "Error in AssignValues: No LesHouches File is found in the given path!"
