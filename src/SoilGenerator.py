import numpy as np
import os
import matplotlib.pyplot as pl
import shutil
from DelPhynMat import DelPhynMat
from vanGenuchten import vanGenuchten
import utils



class SoilGenerator:

    def __init__(self, sand: float, silt: float, clay: float, density: float, van_genuchten: vanGenuchten):
        self.f_sand = sand # m3/m3
        self.f_silt = silt # m3/m3
        self.f_clay = clay # m3/m3
        self.density = density  # kg/m3
        self.soil_type = ""
        self.van_genuchten = van_genuchten


    def __create_mrc_function(self):
        self.h_min = 1e-5  # in m
        self.h_max = 1e12  # in m
        self.h_mvg_trad = np.geomspace(self.h_min, self.h_max, num=10000)
        # --> create function and use fitted polynom for theta/theta_s < transition
        # it might be necessary to modify max_slope_change
        transition = 0.95
        self.theta_eff = self.van_genuchten.theta_s
        self.pc_mvg_trad = 1000 * 9.81 * self.h_mvg_trad  # in Pa
        self.theta_mvg_with_poly = traditional_mvg_with_polynom(self.h_mvg_trad, self.van_genuchten.n, self.van_genuchten.theta_r,
                                                           self.van_genuchten.theta_s, self.van_genuchten.alpha,
                                                           transition=transition)


    def create_m6_file(self, m6_file_path, show_plot=False):
        """
        :param m6_file_path: file path of m6 file
        :param show_plot: if True, the functions will be plotted
        """
        # to SI values
        self.van_genuchten.K0 = self.van_genuchten.K0 / (100 * (24 * 3600))  # K0 [cm/d] -> [m/s]
        self.van_genuchten.alpha = self.van_genuchten.alpha * 100  # alpha [1/cm] -> [1/m]
        porosity = 1 - self.density / 2650

        # determine soil type using DIN table
        self.soil_type = utils.read_soil_type_from_table(self.f_sand, self.f_silt, self.f_clay)

        # Thermal conductivity (model Markert et al.)
        theta_lambda = np.linspace(1e-6, self.van_genuchten.theta_s, num=50)
        lambda_unfrozen, lambda_frozen = thermalConductivity_Markert_modified(theta_lambda, self.f_clay, self.f_sand, self.f_silt,
                                                                                 self.density / 1000, group='all')

        # 1. Water retention curve according to traditional mvg
        self.__create_mrc_function()

        # create sparse spline
        max_slope_change = 1e-4
        idx_sparse = sparse_spline(np.log10(self.h_mvg_trad), self.theta_mvg_with_poly, max_slope_change=max_slope_change)
        theta_mvg_sparse = self.theta_mvg_with_poly[idx_sparse]
        pc_mvg_sparse = self.pc_mvg_trad[idx_sparse]

        # inverse  functions
        theta_mvg_inv_sparse = theta_mvg_sparse[::-1]
        pc_mvg_inv_sparse = pc_mvg_sparse[::-1]

        # 2. Liquid conductivity according to modified MVG, Vogel et al. 2001
        h_Kl = np.geomspace(self.h_min, self.h_max, num=100)  # in m
        theta_Kl, Kl = modified_mvg_Kl(h_Kl, self.van_genuchten.n, self.van_genuchten.theta_r, self.van_genuchten.theta_s,
                                          self.van_genuchten.alpha, self.van_genuchten.K0, self.van_genuchten.l)
        # -> fit to Kleff
        try:
            Kl_eff = KlEff_wessolek(self.soil_type, self.density / 1000, self.f_clay, self.f_silt)
            Kl = Kl * Kl_eff / Kl[0]
        except:
            print("Unknown soil type '{}'.\nCould not calculate correction value for liquid conductivity Kl.".format(self.soil_type))
        # -> clip values which are lgKl < -20
        cut_Kl = Kl < 1e-20
        Kl = Kl[~cut_Kl]
        theta_Kl = theta_Kl[~cut_Kl]

        # plot
        material_name = "{} {} kg/m³".format(self.soil_type, self.density)
        if show_plot:
            fig, ax = pl.subplots(1, 3, figsize=(16, 8))
            ax[0].set_title("Feuchtespeicherfunktion")
            ax[0].semilogx(pc_mvg_sparse, theta_mvg_sparse, 'o-', label=material_name, markersize=2)
            ax[0].set_xlabel("pc [Pa]")
            ax[0].set_ylabel("Theta [m³/m³]")
            ax[0].legend()

            ax[1].set_title("Flüssigwasserleitfähigkeit")
            ax[1].semilogy(theta_Kl, Kl)
            ax[1].set_xlabel("Theta [m³/m³]")
            ax[1].set_ylabel("Kl [s]")

            ax[2].set_title("Wärmeleitfähigkeit")
            ax[2].plot(theta_lambda, lambda_unfrozen, label="unfrozen")
            ax[2].plot(theta_lambda, lambda_frozen, '--', label="frozen")
            ax[2].set_xlabel("Theta [m³/m³]")
            ax[2].set_ylabel("Lambda [W/mK]")
            ax[2].legend()

        # write m6
        current_dir = os.path.dirname(os.path.realpath(__file__))
        m6_dummy_filepath = os.path.join(current_dir, '../data/dummy.m6')
        shutil.copyfile(m6_dummy_filepath, m6_file_path)
        dm = DelPhynMat(m6_file_path)
        dm.writeBaseProperty("NAME", 'DE: ' + material_name)
        dm.writeBaseProperty("RHO", self.density)
        dm.writeBaseProperty("CE", dryheatCapacity_massSpecific())
        dm.writeBaseProperty("MEW", mew_value_fixed())
        dm.writeBaseProperty("THETA_POR", porosity)
        dm.writeBaseProperty("THETA_EFF", self.theta_eff)
        dm.writeBaseProperty("LAMBDA", lambda_unfrozen[0])
        dm.writeBaseProperty("KLEFF", Kl[0])
        dm.writeMaterialFunction(np.log10(pc_mvg_sparse), theta_mvg_sparse, "Theta_l(pC)")
        dm.writeMaterialFunction(theta_Kl[::-1], np.log10(Kl[::-1]), "lgKl(Theta_l)")
        dm.writeMaterialFunction(theta_lambda, lambda_unfrozen, "lambda(Theta_l)")
        dm.writeMaterialFunction(theta_lambda, lambda_frozen, "lambda(Theta_ice)")
        dm.writeMaterialFunction(theta_mvg_inv_sparse, np.log10(pc_mvg_inv_sparse), "pC(Theta_l)")
        # color
        try:
            dm.writeBaseProperty("COLOUR", hex_color(self.soil_type))
        except:
            pass

        if show_plot:
            pl.show()



def thermalConductivity_Markert_modified(theta, fClay, fSand, fSilt, rho, group='all'):
    """
    according to Markert et al. 2017, modified to also include frozen values
    :param theta: in m³/m³
    :param fClay: share of Ton (Clay) in -
    :param fSand: share of Sand in -
    :param fSilt: share of Schluff (Silt) in -
    :param rho: Trockenrohdichte in g/cm3
    :return:
    """
    # parameters from Table 2 of the paper
    if group == 'all':
        p = [1.21, -1.55, 0.02, 0.25, 2.29, 2.12, -1.04, -2.03]
    elif (group is None) and (fSilt + 2* fClay < 0.3) or group=="sand":
        print('group sand chosen')
        p = [1.02, -1.64, -0.13, 0.25, 1.99, 2.87, -1.32, -2.39]
    elif (group is None) and (fSilt > 0.5 and fClay < 0.27) or group == "silt":
        print('group silt chosen')
        p = [1.48, -2.15, 0.78, 0.23, 0.00, 0.86, 0.41, 0.20]
    elif group is None or group=="loam":
        print('group Loam chosen')
        p = [1.64, -2.39, -0.42, 0.28, 3.88, 1.62, -1.10, -2.36]
    else:
        raise Exception('Invalid group')

    rho_si = rho * 1000 # [g/cm3] -> [kg/m3]
    porosity = 1 - rho_si/2650
    lambda_dry = p[0] + p[1] * porosity
    alpha = p[2] * fClay + p[3]
    beta = p[4] * fSand + p[5] * rho + p[6] * fSand * rho + p[7]

    lambda_u = lambda_dry + np.exp(beta - theta ** (-alpha))

    lambda_sat_u = lambda_u[-1]
    lambda_ice = 2.2
    lambda_water = 0.57
    lambda_solid = (lambda_sat_u/lambda_water**porosity)**(1/(1-porosity))
    lambda_sat_f = lambda_solid**(1-porosity) * lambda_ice**porosity
    beta_f = np.log(lambda_sat_f-lambda_dry) + theta[-1]**-alpha
    lambda_f = lambda_dry + np.exp(beta_f - theta ** (-alpha))

    return lambda_u, lambda_f



def dryheatCapacity_massSpecific():
    """
    according to Wang et al 2019
    :return: in J/kgK
    """
    return 750.0


def mew_value_fixed():
    """
    according to IBK measurements and Jabro 2009
    :return: in -
    """
    return 7.7


def traditional_mvg(h, n, theta_r, theta_s, alpha):
    """
    traditional MVG model, water content
    :param h: in m
    :param n: van-Genuchten Parameter in -
    :param theta_r: van-Genuchten Parameter in m3/m3
    :param theta_s: van-Genuchten Parameter in m3/m3
    :param alpha: van-Genuchten Parameter in 1/m
    :return:
    """
    m = 1 - 1/n
    theta_mvg = theta_r + (theta_s - theta_r) / (1 + (alpha * np.abs(h)) ** n) ** m
    return theta_mvg


def traditional_mvg_with_polynom(h, n, theta_r, theta_s, alpha, transition, mod_factor=1.0):
    """
    traditional MVG model, water content
    :param h: in m
    :param n: van-Genuchten Parameter in -
    :param theta_r: van-Genuchten Parameter in m3/m3
    :param theta_s: van-Genuchten Parameter in m3/m3
    :param alpha: van-Genuchten Parameter in 1/m
    :return:
    """
    theta_mvg = traditional_mvg(h, n, theta_r, theta_s, alpha)

    # we fit a 2nd order polynom to the first part to ensure a more steep slope
    if transition >= 1:
        raise Exception('transition must be below 1')

    theta_eff_mod = mod_factor * theta_s
    h_min = np.min(h)
    h_internal = np.geomspace(h_min, np.max(h), num=100000)
    theta_mvg_internal = traditional_mvg(h_internal, n, theta_r, theta_s, alpha)
    theta_trans = transition * theta_s
    slope_mvg = np.gradient(theta_mvg_internal, h_internal)
    h_trans = np.interp(theta_trans, theta_mvg_internal[::-1], h_internal[::-1])
    slope_theta_trans = np.interp(h_trans, h_internal, slope_mvg) # slope at transition point
    # setup LES and solve
    A = [[h_min ** 2, h_min, 1], [h_trans ** 2, h_trans, 1], [2 * h_trans, 1, 0]]
    b = [theta_eff_mod, theta_trans, slope_theta_trans]
    p = np.linalg.solve(A, b)
    # polynom
    theta_poly = p[0] * h**2 + p[1] * h + p[2]
    # combined function
    theta_mvg[h<h_trans] = theta_poly[h<h_trans]

    return theta_mvg


def modified_mvg_Kl(h, n, theta_r, theta_s, alpha, K0, l):
    """
     modified MVG, Vogel et al. 2001, water conductivity
    :param h: in m
    :param n: van-Genuchten Parameter in -
    :param theta_r: van-Genuchten Parameter in m3/m3
    :param theta_s: van-Genuchten Parameter in m3/m3
    :param alpha: van-Genuchten Parameter in 1/m
    :param K0: van-Genuchten Parameter in m/s
    :param l: van-Genuchten Parameter in -
    :return:
    """
    hs = 0.04 # in m
    m = 1 - 1 / n
    theta_m = theta_r + (theta_s - theta_r) * (1 + (alpha * np.abs(hs))**n )**m
    theta_Kl = theta_r + (theta_m - theta_r) / (1 + (alpha * np.abs(h)) ** n) ** m
    theta_Kl[h<=hs] = theta_s
    Se3 = (theta_Kl - theta_r) / (theta_s - theta_r)
    Se_star = lambda Se: (theta_s - theta_r) / (theta_m - theta_r) * Se
    F = lambda Se_star: (1 - Se_star ** (1 / m)) ** m
    Km = K0 * Se3 ** l * ((1 - F(Se_star(Se3))) / (1 - F(Se_star(1)) ))**2
    Km[h<=hs] = K0
    Kl = Km / 9.81

    return theta_Kl, Kl


def KlEff_wessolek(soil_type, rho, fClay, fSilt):
    """
    according to wessolek rote reihe 2, tabelle 8
    :param fSilt: in -
    :param fClay: in -
    :param soil_type:
    :param rho: in g/cm³
    :return:
    """

    # Lagerungsdichte according to rote reihe 1, Seite 6
    ld = rho + 0.005 * fClay*100 + 0.001 * fSilt*100
    print('Lagerungsdichte = {:.2f}'.format(ld))

    # according to wessolek rote reihe 2, tabelle 8
    if soil_type in ["Ss", "mS"]:
        Kf = 1384 - 974 * ld + 159 * ld ** 2
    elif soil_type == "fS":
        Kf = 1452 - 1330 * ld + 1312 * ld ** 2
    elif soil_type == "gS":
        Kf = 7840 - 7895 * ld + 2030.5 * ld ** 2
    elif soil_type == "Uu":
        Kf = 207 - 184 * ld + 40.9 * ld ** 2
    elif soil_type == "Ut2":
        Kf = 250 - 220 * ld + 48 * ld ** 2
    elif soil_type == "Ut3":
        Kf = 262 - 237.7 * ld + 52.8 * ld ** 2
    elif soil_type == "Ut4":
        Kf = 288.7 - 271.5 * ld + 63.9 * ld ** 2
    elif soil_type == "Us":
        Kf = 240.9 - 228.1 * ld + 54.2 * ld ** 2
    elif soil_type == "Uls":
        Kf = 199 - 158 * ld + 29.6 * ld ** 2
    elif soil_type == "Lu":
        Kf = 218.5 - 181 * ld + 37 * ld ** 2
    elif soil_type in ["Tt", "Tu2", "Tl", "Tu3"]:
        Kf = 527 - 535 * ld + 135.9 * ld**2
    elif soil_type == "Tu4":
        Kf = 501 - 499.7 * ld + 125 * ld ** 2
    elif soil_type == "Lt3":
        Kf = 258.3 - 222 * ld + 47.1 * ld ** 2
    elif soil_type == "Lt2":
        Kf = 262 - 184.9 * ld + 27.8 * ld ** 2
    elif soil_type == "Lts":
        Kf = 218.4 - 181 * ld + 37 * ld ** 2
    elif soil_type in ["Ts2", "Ts3"]:
        Kf = 390.4 - 358.3 * ld + 82.7 * ld ** 2
    elif soil_type == "Ts4":
        Kf = 387.5 - 331.4 * ld + 70.4 * ld ** 2
    elif soil_type == "Sl2":
        Kf = 1022 - 856 * ld + 180.3 * ld**2
    elif soil_type == "Sl3":
        Kf = 548 - 471 * ld + 103.2 * ld ** 2
    elif soil_type == "Sl4":
        Kf = 476.5 - 405.7 * ld + 87.4 * ld ** 2
    elif soil_type == "Slu":
        Kf = 383 - 306 * ld + 58.7 * ld ** 2
    elif soil_type == "St2":
        Kf = 765.6 - 596.1 * ld + 114.3 * ld ** 2
    elif soil_type == "St3":
        Kf = 691.2 - 595.8 * ld + 128.9 * ld ** 2
    elif soil_type == "Su2":
        Kf = 995.7 - 895.9 * ld + 205.2 * ld ** 2
    elif soil_type == "Su3":
        Kf = 530 - 454.4 * ld + 97.7 * ld ** 2
    elif soil_type == "Su4":
        Kf = 425.8 - 345.7 * ld + 69 * ld ** 2
    elif soil_type == "Ls2":
        Kf = 265 - 201.6 * ld + 36 * ld ** 2
    elif soil_type == "Ls3":
        Kf = 310.2 - 248 * ld + 49 * ld ** 2
    elif soil_type == "Ls4":
        Kf = 240 - 133.7 * ld + 9.5 * ld ** 2
    else:
        raise Exception("Unknown soil type!")

    # cm/d -> m/s
    Kf_si = Kf / (100 * (24*3600))
    Kl_eff = Kf_si/9.81

    return Kl_eff


def sparse_spline(x, y, max_slope_change=0.1):
    # consider the slope of x,y: where slope changes more than max_slope_change, index points are generated, returns the index points
    slope = np.gradient(y, x)
    idx_sparse = [0,]
    prev_slope = slope[0]
    for n, s in enumerate(slope[1:]):
        slope_change = (prev_slope - slope[n])
        if np.abs(slope_change) > max_slope_change:
            idx_sparse.append(n)
            prev_slope = slope[n]

    if len(x)-1 not in idx_sparse:
        idx_sparse.append(len(x)-1)

    return idx_sparse


def hex_color(soil_type):
    if soil_type in ['fS', 'mS', 'Ss', 'gS']:
        return '#fffdd6'
    elif soil_type in ['St2', 'Su2', 'Sl2', 'Sl3']:
        return '#feff8a'
    elif soil_type in ['Su3', 'Su4']:
        return '#fef600'
    elif soil_type in ['Slu', 'Sl4', 'St3']:
        return '#da9e60'
    elif soil_type in ['Lt2', 'Ls2', 'Ls3', 'Ls4']:
        return '#b86942'
    elif soil_type in ['Lts', 'Ts3', 'Ts4']:
        return '#bb8141'
    elif soil_type in ['Us', 'Uu']:
        return '#f1c5ac'
    elif soil_type in ['Ut2', 'Ut3', 'Uls']:
        return '#eea481'
    elif soil_type in ['Ut4', 'Lu']:
        return '#ea7846'
    elif soil_type in ['Tu3', 'Tu4', 'Lt3']:
        return '#f39dc0'
    elif soil_type in ['Tt', 'Tu2', 'Tl', 'Ts2']:
        return '#994775'
    else:
        raise Exception('soil type \'{}\' not found'.format(soil_type))



