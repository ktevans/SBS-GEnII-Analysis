import ROOT
from ROOT import RDataFrame
import sys
import math
import uproot
import numpy as np
#from Scientific.Geometry import Vector
import pandas as pd
import matplotlib.pyplot as plt

## First we get the file names sorted
ConfigName = sys.argv[1]
#/cache/halla/sbs/prod/GEnII/pass1/GEN2/He3/rootfiles/
#GEN2/He3/rootfiles/e1209016_fullreplay_2321_stream0_2_seg3_3
#.root

inFileName = "/cache/halla/sbs/prod/GEnII/pass1/" + str(ConfigName) + ".root"

kin = ""
target = ""
Q2_kin = 0.0 #Q2 expected from kinematic in GeV2
p_electron = 0.0 #electron momentum in GeV
p_nucleon = 0.0 #nucleon momentum in GeV
energy_beam = 0.0 #beam energy expected from kinematic
bb_angle = 0.0 #bigbite magnet angle (what angle?) in degrees
bb_dist = 1.63 #distance from center of target to bb magnet in meters
sbs_angle = 0.0 #super bigbite magnet angle (what angle?) in degrees
sbs_dist = 2.8 #distance from center of target to sbs magnet in meters
hcal_angle = 0.0 #hcal angle (what angle?) in degrees
hcal_dist = 17.0 #distance from center of target to hcal in meters

if "GEN2" in ConfigName:
    kin = "GEN2"
    Q2_kin = 3.00
    p_electron = 2.69
    p_nucleon = 2.37
    energy_beam = 4.291
    bb_angle = 29.5
    sbs_angle = 34.7
    hcal_angle = 34.7
if "GEN3" in ConfigName:
    kin = "GEN3"
    Q2_kin = 6.83
    p_electron = 2.73
    p_nucleon = 4.51
    energy_beam = 6.373
    bb_angle = 36.5
    sbs_angle = 22.1
    hcal_angle = 21.6
if "GEN4a" in ConfigName:
    kin = "GEN4a"
    Q2_kin = 9.82
    p_electron = 3.21
    p_nucleon = 6.11
    energy_beam = 8.448
    bb_angle = 35.0
    sbs_angle = 18.0
    hcal_angle = 18.0
if "GEN4b" in ConfigName:
    kin = "GEN4b"
    Q2_kin = 9.82
    p_electron = 3.21
    p_nucleon = 6.11
    energy_beam = 8.448
    bb_angle = 35.0
    sbs_angle = 18.0
    hcal_angle = 18.0

if "He3" in ConfigName:
    target = "He3"
if "H2" in ConfigName:
    target = "H2"

outFileName = "/volatile/halla/sbs/ktevans/" + str(kin) + "_" + str(target) + "_dxdy.root"

print("Reading from " + str(inFileName) + " and writing to " + str(outFileName))

## Physics constants:

PI = np.pi
c = 299792458 # m/s -> speed of light
q_e = 0.0000000000000000001602176634 #1.602176634e^-19 C -> charge of electron
M_e = 0.000510998950 # GeV/c^2 -> mass of electron
M_p = 0.93827208816 # GeV/c^2 -> mass of proton
M_n = 0.93956542052 # GeV/c^2 -> mass of neutron
alpha = 0.00729927 # unitless -> fine structure constant

## Next we define some methods that will help us.

def ListBranches(fileName):
    ## This method will list out all the branch names in a root file. It's typically useful initially, but not needed often.

    inFile = ROOT.TFile(fileName)
    t = inFile.T

    for branch in t.GetListOfBranches():
        if "Ndata" not in branch.GetName():
            print("\"%s\", " %branch.GetName(), end="")

class dxdy:

    def __init__(self):
        self.treeData = pd.DataFrame()

    def GetNumpyArray(self,fileName):
    ## This method will make an array out of the TTree of a standard SBS root file.
        file = uproot.open(fileName)
        T = file["T"]

        data = T.arrays(["bb.bbtrig.a_amp_p", "bb.bbtrig.a_time", "bb.bbtrig.adcelemID", "bb.bbtrig.tdc", "bb.bbtrig.tdcelemID", "bb.eps_over_etot", "bb.etot_over_p","bb.grinch_tdc.allclus.adc", "bb.grinch_tdc.allclus.size", "bb.grinch_tdc.allclus.t_mean", "bb.grinch_tdc.allclus.t_rms", "bb.grinch_tdc.allclus.tot_mean", "bb.grinch_tdc.allclus.trackindex", "bb.grinch_tdc.allclus.x_mean", "bb.grinch_tdc.allclus.y_mean", "bb.grinch_tdc.hit.amp", "bb.grinch_tdc.hit.clustindex", "bb.grinch_tdc.hit.col", "bb.grinch_tdc.hit.pmtnum", "bb.grinch_tdc.hit.row", "bb.grinch_tdc.hit.time", "bb.grinch_tdc.hit.trackindex", "bb.grinch_tdc.hit.xhit", "bb.grinch_tdc.hit.yhit", "bb.hodotdc.clus.bar.tdc.id", "bb.hodotdc.clus.bar.tdc.itrack", "bb.hodotdc.clus.bar.tdc.meantime", "bb.hodotdc.clus.bar.tdc.meantot", "bb.hodotdc.clus.bar.tdc.timediff", "bb.hodotdc.clus.bar.tdc.timehitpos", "bb.hodotdc.clus.bar.tdc.vpos", "bb.hodotdc.clus.id", "bb.hodotdc.clus.size", "bb.hodotdc.clus.tdiff", "bb.hodotdc.clus.tmean", "bb.hodotdc.clus.totmean", "bb.hodotdc.clus.trackindex", "bb.hodotdc.clus.xmean", "bb.hodotdc.clus.ymean", "bb.ps.clus.again", "bb.ps.clus.atime", "bb.ps.clus.col", "bb.ps.clus.e", "bb.ps.clus.eblk", "bb.ps.clus.id", "bb.ps.clus.nblk", "bb.ps.clus.row", "bb.ps.clus.tdctime", "bb.ps.clus.x", "bb.ps.clus.y", "bb.ps.clus_blk.again", "bb.ps.clus_blk.atime", "bb.ps.clus_blk.col", "bb.ps.clus_blk.e", "bb.ps.clus_blk.id", "bb.ps.clus_blk.row", "bb.ps.clus_blk.tdctime", "bb.ps.clus_blk.x", "bb.ps.clus_blk.y", "bb.sh.clus.again", "bb.sh.clus.atime", "bb.sh.clus.col", "bb.sh.clus.e", "bb.sh.clus.eblk", "bb.sh.clus.id", "bb.sh.clus.nblk", "bb.sh.clus.row", "bb.sh.clus.tdctime", "bb.sh.clus.x", "bb.sh.clus.y", "bb.sh.clus_blk.again", "bb.sh.clus_blk.atime", "bb.sh.clus_blk.col", "bb.sh.clus_blk.e", "bb.sh.clus_blk.id", "bb.sh.clus_blk.row", "bb.sh.clus_blk.tdctime", "bb.sh.clus_blk.x", "bb.sh.clus_blk.y", "bb.tdctrig.tdc", "bb.tdctrig.tdcelemID", "bb.tr.beta", "bb.tr.chi2", "bb.tr.d_ph", "bb.tr.d_th", "bb.tr.d_x", "bb.tr.d_y", "bb.tr.dbeta", "bb.tr.dtime", "bb.tr.flag", "bb.tr.ndof", "bb.tr.p", "bb.tr.pathl", "bb.tr.ph", "bb.tr.px", "bb.tr.py", "bb.tr.pz", "bb.tr.r_ph", "bb.tr.r_th", "bb.tr.r_x", "bb.tr.r_y", "bb.tr.tg_dp", "bb.tr.tg_ph", "bb.tr.tg_th", "bb.tr.tg_x", "bb.tr.tg_y", "bb.tr.th", "bb.tr.time", "bb.tr.vx", "bb.tr.vy", "bb.tr.vz", "bb.tr.x", "bb.tr.y", "bb.x_bcp", "bb.x_fcp", "bb.y_bcp", "bb.y_fcp", "bb.z_bcp", "bb.z_fcp", "sbs.hcal.clus.again", "sbs.hcal.clus.atime", "sbs.hcal.clus.col", "sbs.hcal.clus.e", "sbs.hcal.clus.eblk", "sbs.hcal.clus.id", "sbs.hcal.clus.nblk", "sbs.hcal.clus.row", "sbs.hcal.clus.tdctime", "sbs.hcal.clus.x", "sbs.hcal.clus.y", "sbs.hcal.clus_blk.again", "sbs.hcal.clus_blk.atime", "sbs.hcal.clus_blk.col", "sbs.hcal.clus_blk.e", "sbs.hcal.clus_blk.id", "sbs.hcal.clus_blk.row", "sbs.hcal.clus_blk.tdctime", "sbs.hcal.clus_blk.x", "sbs.hcal.clus_blk.y", "sbs.trig.a_amp_p", "sbs.trig.a_p", "sbs.trig.a_time", "sbs.trig.adcelemID","bb.grinch_tdc.bestcluster", "bb.grinch_tdc.clus.adc", "bb.grinch_tdc.clus.size", "bb.grinch_tdc.clus.t_mean", "bb.grinch_tdc.clus.t_rms", "bb.grinch_tdc.clus.tot_mean", "bb.grinch_tdc.clus.trackindex", "bb.grinch_tdc.clus.x_mean", "bb.grinch_tdc.clus.y_mean", "bb.grinch_tdc.nclus", "bb.grinch_tdc.ngoodhits", "bb.grinch_tdc.ntrackmatch", "bb.hodotdc.nclus", "bb.ps.againblk", "bb.ps.atimeblk", "bb.ps.colblk", "bb.ps.e", "bb.ps.eblk", "bb.ps.idblk", "bb.ps.index", "bb.ps.nblk", "bb.ps.nclus", "bb.ps.ngoodADChits", "bb.ps.rowblk", "bb.ps.x", "bb.ps.y", "bb.sh.againblk", "bb.sh.atimeblk", "bb.sh.colblk", "bb.sh.e", "bb.sh.eblk", "bb.sh.idblk", "bb.sh.index", "bb.sh.nblk", "bb.sh.nclus", "bb.sh.ngoodADChits", "bb.sh.rowblk", "bb.sh.x", "bb.sh.y", "bb.tr.n", "bb.ts.over_threshold", "e.kine.Q2", "e.kine.W2", "e.kine.angle", "e.kine.epsilon", "e.kine.nu", "e.kine.omega", "e.kine.ph_q", "e.kine.q3m", "e.kine.q_x", "e.kine.q_y", "e.kine.q_z", "e.kine.th_q", "e.kine.x_bj", "g.datatype", "g.evlen", "g.evnum", "g.evtime", "g.evtyp", "g.runnum", "g.runtime", "g.runtype", "g.trigbits", "sbs.hcal.againblk", "sbs.hcal.atimeblk", "sbs.hcal.colblk", "sbs.hcal.e", "sbs.hcal.eblk", "sbs.hcal.idblk", "sbs.hcal.index", "sbs.hcal.nblk", "sbs.hcal.nclus", "sbs.hcal.ngoodADChits", "sbs.hcal.rowblk", "sbs.hcal.tdctimeblk", "sbs.hcal.x", "sbs.hcal.y"], library="pd")

        self.treeData = data

    def DrawTestHisto(self):
        df = pd.DataFrame()
        df = self.treeData

        plt.scatter(df["bb.ps.e"],df["bb.sh.e"])
        plt.show()

    def ExpectedValuesMethod3(self):
        df = pd.DataFrame()
        df = self.treeData

        hcal_xaxis = np.array([0, -1, 0])
        hcal_zaxis = np.array([np.sin(np.radians(-hcal_angle)), 0, np.cos(np.radians(-hcal_angle))])
        hcal_yaxis = np.cross(hcal_zaxis, hcal_xaxis)

        hcal_offset = 0.1 #change later
        hcal_origin = (hcal_dist * hcal_zaxis) + (hcal_offset * hcal_xaxis)

        #print(hcal_xaxis)
        #print(hcal_zaxis)
        #print(hcal_yaxis)
        #print(hcal_origin)

        df["bb.tr.vz"] = df["bb.tr.vz"].apply(lambda x: x[0])
        df["vertex_0"] = 0.0
        df["vertex_1"] = 0.0
        df["vertex"] = df[["vertex_0", "vertex_1", "bb.tr.vz"]].values.tolist()

        df.drop("vertex_0", axis=1, inplace=True)
        df.drop("vertex_1", axis=1, inplace=True)

        p_beam_array = np.array([0, 0, energy_beam, energy_beam])
        df["p_beam"] = [p_beam_array] * len(df)

        df["bb.tr.px"] = df["bb.tr.px"].apply(lambda x: x[0])
        df["bb.tr.py"] = df["bb.tr.py"].apply(lambda x: x[0])
        df["bb.tr.pz"] = df["bb.tr.pz"].apply(lambda x: x[0])
        df["bb.tr.p"] = df["bb.tr.p"].apply(lambda x: x[0])

        df["p_eprime"] = df[["bb.tr.px","bb.tr.py","bb.tr.pz","bb.tr.p"]].values.tolist()
        
        p_targ_array = np.array([0, 0, 0, 0.5*(M_n+M_p)])
        df["p_targ"] = [p_targ_array] * len(df)

        df["q"] = df["p_beam"] - df["p_eprime"]

        df["e_theta"] = df["bb.tr.pz"] / df["bb.tr.p"]

        Nmass = 0.5 * (M_p + M_n)
        df["p_central"] = energy_beam / (1.0 + (energy_beam / Nmass) * (1.0 - np.cos(df["e_theta"])))


if __name__ == '__main__':

    calculation = dxdy()
    calculation.GetNumpyArray(str(inFileName))
#    calculation.DrawTestHisto()
    calculation.ExpectedValuesMethod3()
