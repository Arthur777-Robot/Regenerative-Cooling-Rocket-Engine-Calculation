# -*- coding: utf-8 -*-
'''再生冷却設計ツール
'''
import sys
reload(sys)
import platform
# デフォルトの文字コードを変更する．
sys.setdefaultencoding('utf-8')
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib.font_manager import FontProperties
import xlrd
import csv
import pandas as pd
import pprint

if 'Windows' == platform.system():
    font_prop = FontProperties(fname=r'C:\WINDOWS\Fonts\MSGothic.ttf')

if 'Darwin' == platform.system(): # for Mac
    font_prop = FontProperties(fname='/Library/Fonts/ipaexg.ttf')
    matplotlib.rcParams['font.family'] = font_prop.get_name()
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['savefig.dpi'] = 300
    matplotlib.rcParams['mathtext.default'] = 'regular'

plt.close('all')

class Parameter():
    def __init__(self, file_name = u"再生冷却.xlsx"):
        print('Paramter...')
        book = file_name
        sheet = "Input"
        EXL = pd.ExcelFile(book)
        self.df = EXL.parse(sheet)
        self.thrust = self.df.data[u"推力"]
        self.press_chamber = self.df.data[u"燃焼室圧力"]
        self.of_ratio = self.df.data[u"O/F"]
        self.type_fuel = self.df.data[u"推進剤（燃料）"]
        self.type_oxidizer = self.df.data[u"推進剤（酸化剤）"]
        self.density_fuel = self.df.data[u"燃料密度"]
        self.density_oxidizer = self.df.data[u"酸化剤密度"]
        self.L_star = self.df.data[u"L*"]
        self.tensile_strength = self.df.data[u"燃焼室強度"]
        self.safety_coef = self.df.data[u"安全率"]
        self.thermal_capacity_fuel = self.df.data[u"燃料熱容量"]
        self.visc_fuel = self.df.data[u"燃料粘性係数"]
        self.thermal_conductivity_fuel = self.df.data[u"燃料熱伝導率"]
        self.thermal_conductivity_metal = self.df.data[u"金属熱伝導率"]

    def __repr__(self):
        pp = pprint.PrettyPrinter(indent=4)
        return str(pp.pprint(self.__dict__))

    def run_CEA(self):
        file_name = '00temp' # 一時的に作られるファイル
        # ---- .inpファイル作り ----
        input_file = file_name + '.inp'
        str = """problem    o/f=%.1f,
            rocket  equilibrium  frozen  nfz=2
          p,bar = %d
          pi/p= %d
        react
          fuel=%s  wt=100  t,k=298.15
          oxid=%s wt=100  t,k=90.17
        output transport
        output short
        output trace=1e-5
            plot p pi/p o/f isp ivac aeat cf t
        end
        """ % (self.of_ratio, self.press_chamber*10,
               self.press_chamber*10,
               self.type_fuel, self.type_oxidizer)
        f = open(input_file, 'w') # 書き込みモードで開く
        f.write(str) # 引数の文字列をファイルに書き込む
        f.close() # ファイルを閉じる
        # ---- FCEA2コマンド ----
        cmd = 'FCEA2' # FCEA2にはパスを通しておく
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        p.communicate(file_name)

    def read_CEA(self):
        col_name = ['c{0:02d}'.format(i) for i in range(20)]
        df = pd.read_table("00temp.out", skiprows=10,
                           names = col_name, delim_whitespace=True)
        # データ読み込みする
        # np.array("CHAMBER", "THROAT", "EXIT")
        # _e:equilibrium, _f:frozen
        date_name = ["T,","VISC,MILLIPOISE","Cp,","PRANDTL","Ae/At","CSTAR,","CF","Isp,"]
        self.temp_chamber_e = np.array(df[df.c00 == date_name[0]][0:1][['c02','c03','c04']], dtype=np.float64)
        self.temp_chamber_f = np.array(df[df.c00 == date_name[0]][1:2][['c02','c03','c04']], dtype=np.float64)
        self.visc_gas_e = np.array(df[df.c00 == date_name[1]][0:1][['c01','c02','c03']], dtype=np.float64)
        self.visc_gas_f = np.array(df[df.c00 == date_name[1]][1:2][['c01','c02','c03']], dtype=np.float64)
        self.Cp_e = np.array(df[df.c00 == date_name[2]][0:1][['c02','c03','c04']], dtype=np.float64)
        self.Cp_f = np.array(df[df.c00 == date_name[2]][3:4][['c02','c03','c04']], dtype=np.float64)
        self.prandtl_e = np.array(df[df.c00 == date_name[3]][0:1][['c02','c03','c04']], dtype=np.float64)
        self.prandtl_f = np.array(df[df.c00 == date_name[3]][2:3][['c02','c03','c04']], dtype=np.float64)
        self.throat_ratio_e = np.array(df[df.c00 == date_name[4]][0:1][['c01','c02']], dtype=np.float64)[0][1]
        self.throat_ratio_f = np.array(df[df.c00 == date_name[4]][1:2][['c01','c02']], dtype=np.float64)[0][1]
        self.C_star_e = np.array(df[df.c00 == date_name[5]][0:1][['c02','c03']], dtype=np.float64)[0][1]
        self.C_star_f = np.array(df[df.c00 == date_name[5]][1:2][['c02','c03']], dtype=np.float64)[0][1]
        self.Cf_e = np.array(df[df.c00 == date_name[6]][0:1][['c01','c02']], dtype=np.float64)[0][1]
        self.Cf_f = np.array(df[df.c00 == date_name[6]][1:2][['c01','c02']], dtype=np.float64)[0][1]
        self.Isp_e = np.array(df[df.c00 == date_name[7]][0:1][['c02','c03']], dtype=np.float64)[0][1]
        self.Isp_f = np.array(df[df.c00 == date_name[7]][1:2][['c02','c03']], dtype=np.float64)[0][1]

    def select_equilibrium_or_frozen(self):
        print("calculate with frozen or equilibrium option\n if frozen, enter f. if equilibrium, enter e")
        input_char = raw_input()
        if(input_char == 'e'):
            print('/// select equilibrium ///')
            self.temp_chamber = self.temp_chamber_e
            self.visc_gas = self.visc_gas_e
            self.Cp = self.Cp_e
            self.prandtl = self.prandtl_e
            self.throat_ratio = self.throat_ratio_e
            self.C_star = self.C_star_e
            self.Cf = self.Cf_e
            self.Isp = self.Isp_e
        elif(input_char == 'f'):
            print('/// select frozen ///')
            self.temp_chamber = self.temp_chamber_f
            self.visc_gas = self.visc_gas_f
            self.Cp = self.Cp_f
            self.prandtl = self.prandtl_f
            self.throat_ratio = self.throat_ratio_f
            self.C_star = self.C_star_f
            self.Cf = self.Cf_f
            self.Isp = self.Isp_f
        else:
            print('*** select Error ***')

    def calc_nozzle(self):
        self.At = self.thrust / self.press_chamber * self.Cf
        self.Dt = np.sqrt(4.0 * self.At / np.pi)
        self.Ae = self.At * self.throat_ratio
        self.De = np.sqrt(4.0 * self.Ae / np.pi)

    def calc_chamber(self):
        self.Vc = self.At * self.L_star * 1000
        while(1):
            print("input chamber diameter in mm. when satisfied, input type ok")
            print "Dc = ",
            input_temp = raw_input()
            if(input_temp == "ok"):
                break
            elif(input_temp > 0):
                self.Dc = float(input_temp)
                self.Ac = self.Dc **2 * np.pi /4.0
                self.Lc = self.Vc / self.Ac
                print("Lc = %f[mm]" % (self.Lc))
                print("Ac = %f[mm^2]" % (self.Ac))
            else:
                print("** input Error **")

    def calc_chamber_strength(self):
        self.Ct = (self.press_chamber * self.Dc) / \
                  ((2.0 * self.tensile_strength / self.safety_coef) - \
                  1.2 * self.press_chamber)

    def calc_fuel_consumption(self):
        self.mt = self.thrust / self.Isp
        self.mf = self.mt / (self.of_ratio + 1)
        self.mo = self.mt * self.of_ratio / (self.of_ratio + 1)

    def calc_gas_heat_stansfer_coef(self):
        bartz = np.zeros(5)
        bartz[0] = 0.026 * (self.Dt/1000) ** 0.2
        bartz[1] = ((self.visc_gas[0][0]/10000) ** 0.2) * self.Cp[0][0] * 1000 / \
                   (self.prandtl[0][0] ** 0.6)
        bartz[2] = (self.press_chamber * 10e6 / self.C_star) ** 0.8
        bartz[3] = (self.Dt / (1.5 * self.Dt / 2)) ** 0.1
        bartz[4] = (self.At / self.Ac) ** 0.9
        self.hg = bartz.prod()
        self.rg = 1.0 / self.hg

    def calc_fuel_heat_stansfer_coef(self):
        while(1):
            print("input fuel path dimention for cooling")
            print "path width [mm] = ",
            self.path_width = input()
            print "path height [mm] = ",
            self.path_height = input()
            print "number of channel = ",
            self.path_num = input()
            self.path_area = self.path_width * self.path_height
            self.hydraulic_diam = 2.0 * self.path_area / \
                                  (self.path_width + self.path_height)
            self.vf = self.mf / self.density_fuel / 1000.0 / \
                      (self.path_area / 10e6) / self.path_num
            self.delta_p = 0.5 * self.density_fuel * 1000 * self.vf ** 2
            print("path area = %f[mm^2]" % (self.path_area))
            print("hydraulic diameter = %f[mm]" % (self.hydraulic_diam))
            print("fuel velocity in cooling channel = %f[m/sec]" % (self.vf))
            print("channel pressure drop = %f[MPa]" % (self.delta_p / 10.0e6))
            print("if you are satisfied, type ok")
            input_temp = raw_input()
            if(input_temp == "ok"):
                break
        self.nyuf = self.visc_fuel / self.density_fuel / 1000.0
        self.Re = self.vf * self.hydraulic_diam * 1000 / self.nyuf
        self.Pr = self.visc_fuel * self.thermal_capacity_fuel
        self.Nu = 0.023 * (self.Re ** 0.8) * (self.Pr ** 0.6)
        self.hf = self.Nu * self.thermal_conductivity_fuel / \
                  self.hydraulic_diam / 10e6
        self.rf = 1.0 / self.hf

    def calc_metal_heat_stansfer_coef(self):
        self.hm = self.thermal_conductivity_metal
        self.rm = (self.Ct * 10e-3) / self.hm

    def calc_total_heat_stansfer_coef(self):
        self.rc = 0.000407663;	# carbon heat resistivity
        self.rt = self.rg + self.rf + self.rm + self.rc
        self.Tc = self.temp_chamber[0][0]
        self.Q = (self.Tc - 298) / self.rt
        self.Tcw = self.Tc - self.Q * (self.rg + self.rc)
        self.Tcwc = self.Tc - self.Q * (self.rg + self.rc + self.rm)

    def calc_delta_fuel_temp(self):
        self.A_ct = np.pi * (self.Dc + 2 * self.Ct) * self.Lc
        self.A_ct = self.A_ct * 1.2
        self.Tf_out = (self.A_ct * 10e-6) * self.Q / (self.vf * self.thermal_capacity_fuel)

if __name__ == "__main__":
    print("Regenetive cooling...")
    p = Parameter(u"再生冷却.xlsx")
    p.run_CEA()
    p.read_CEA()
    p.select_equilibrium_or_frozen()
    p.calc_nozzle()
    p.calc_chamber()
    p.calc_chamber_strength()
    p.calc_fuel_consumption()
    p.calc_gas_heat_stansfer_coef()
    p.calc_fuel_heat_stansfer_coef()
    p.calc_metal_heat_stansfer_coef()
    p.calc_total_heat_stansfer_coef()
    p.calc_delta_fuel_temp()
