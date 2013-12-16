import time

startTime = time.time()
timeCheckStart = time.time()
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import numpy as np
from scipy.integrate import odeint
print "Import: ", time.time() - timeCheckStart

#Bogus values
temp = np.array([28.8,28.2,26.0,23.8,23.0,22.2,20.0,16.6,14.1,12.2,11.6,10.6,10.8,10.6,10.6,9.4,7.9,4.0,3.0,3.0,3.6,3.8,3.2,2.8,3.1,3.6,2.0,1.8,0.2,-1.5,-4.0,-4.1,-10.3,-11.3,-13.4,-14.1,-15.7,-23.1,-23.3,-24.2,-24.7,-27.5,-31.3,-31.7,-34.1,-37.3,-37.3,-37.5,-37.8,-45.2,-46.3,-47.3,-47.1,-51.5,-51.9,-52.9,-54.5,-54.3,-53.1,-60.0,-60.5,-60.9,-62.0,-62.9,-63.0,-63.1,-66.4,-68.7,-70.9,-70.3])
dewp = np.array([3.8,-3.8,-3.0,-2.2,-2.5,-2.8,-1.9,-0.4,-0.6,-0.8,-1.4,-2.4,-6.2,-5.4,-5.3,0.4,0.3,0.0,-0.1,-1.0,-9.4,-15.2,-11.8,-15.2,-17.8,-22.4,-21.0,-27.2,-29.7,-32.5,-25.3,-25.1,-37.3,-36.3,-36.9,-37.1,-42.7,-41.1,-41.3,-39.5,-38.7,-48.5,-55.3,-59.7,-63.4,-68.3,-68.3,-69.5,-69.7,-75.5,-76.3,-77.3,-77.1,-75.5,-74.9,-74.9,-75.5,-75.5,-75.1,-77.3,-77.5,-76.9,-77.5,-77.9,-78.0,-78.1,-80.2,-81.7,-82.9,-82.3])
pres = np.array([978.0,974.0,949.0,925.0,916.6,908.0,884.7,850.0,823.6,804.0,794.3,780.0,772.0,766.0,765.9,753.0,738.1,700.0,692.0,688.0,685.0,680.0,672.0,662.0,659.9,656.0,638.0,630.0,611.8,592.0,566.6,566.0,506.0,500.0,484.1,479.0,462.0,402.0,400.0,394.6,392.0,373.0,350.0,332.0,318.4,301.0,300.0,293.0,291.7,254.9,250.0,243.0,237.0,215.0,205.0,200.0,185.0,183.8,176.0,151.7,150.0,149.0,144.4,141.0,137.5,136.0,124.4,117.0,106.0,100.0])
wndU = np.array([28.8,28.2,26.0,23.8,23.0,22.2,20.0,16.6,14.1,12.2,11.6,10.6,10.8,10.6,10.6,9.4,7.9,4.0,3.0,3.0,3.6,3.8,3.2,2.8,3.1,3.6,2.0,1.8,0.2,-1.5,-4.0,-4.1,-10.3,-11.3,-13.4,-14.1,-15.7,-23.1,-23.3,-24.2,-24.7,-27.5,-31.3,-31.7,-34.1,-37.3,-37.3,-37.5,-37.8,-45.2,-46.3,-47.3,-47.1,-51.5,-51.9,-52.9,-54.5,-54.3,-53.1,-60.0,-60.5,-60.9,-62.0,-62.9,-63.0,-63.1,-66.4,-68.7,-70.9,-70.3])
wndV = np.array([3.8,-3.8,-3.0,-2.2,-2.5,-2.8,-1.9,-0.4,-0.6,-0.8,-1.4,-2.4,-6.2,-5.4,-5.3,0.4,0.3,0.0,-0.1,-1.0,-9.4,-15.2,-11.8,-15.2,-17.8,-22.4,-21.0,-27.2,-29.7,-32.5,-25.3,-25.1,-37.3,-36.3,-36.9,-37.1,-42.7,-41.1,-41.3,-39.5,-38.7,-48.5,-55.3,-59.7,-63.4,-68.3,-68.3,-69.5,-69.7,-75.5,-76.3,-77.3,-77.1,-75.5,-74.9,-74.9,-75.5,-75.5,-75.1,-77.3,-77.5,-76.9,-77.5,-77.9,-78.0,-78.1,-80.2,-81.7,-82.9,-82.3])

class skewTLogP():
    #Assumes pressure in mb and temp in C
    #Initialization Function
    def __init__(self):
        timeCheckStart = time.time()

        #USER SPECIFIED VALUES (Mostly - I might change what they can do later)
        self.tMin = -40.
        self.tMax = 55.
        self.pMax = 1100.
        self.pMin = 100.
        
        #Plot title
        self.title = 'Skew-T Log-P Test Image'
        
        #Whether or not to create parcel trace info
        parcelTrace = 'T'
        if parcelTrace == 'T':
            pass

        #figSize = (9,12)
        figSize = (10*3/4, 10)

        self.mixRatStop = 800 #Pressure at which to stop plotting mixing ratio
        self.windOffset = 7   #% indent from right bound for wind axis
        self.windOffsetPos = self.tMax - (self.windOffset/100.)*(self.tMax - self.tMin)

#                      TempProf   DewpProf   Isobar     IsoT<0     IsoT=0     IsoT>0     DAdiabat   MAdiabat   MixRat     Wind Axis  Wind Barbs
        self.colors = ['#ff0000', '#009900', '#000000', '#0000ff', '#0000ff', '#ff0000', '#880000', '#009900', '#ff9900', '#000000', '#404040']
        self.lStyle = ['-'      , '-'      , '-'      , '-'      , '-'      , '-'      , '--'     , '--'     , ':'      , '-'      , '-'      ]
        self.lWidth = [1.8      , 1.8      , 1.       , 0.25     , 0.75     , 0.25     , 0.5      , 0.5      , 2.       , 1.       , 0.5      ]
#                      Red        DGreen     Black      Red        Blue       Blue       Dark Red   DGreen     Orangish   Black      Lght Black
#                      Solid      Solid      Solid      Dashed     Dashed     Dashed     Dashed     Dashed     Dotted     Solid      Solid
#                      Normal     Normal     Normal     ExSmall    Med        ExSmall    Small      Small      Large      Normal     Small


        #DATA VALUES - SHOULD BE PROVIDED BY GET DATA FUNCTION
        self.data=np.zeros((len(pres),), dtype=[('pres','f4'), ('temp','f4'), ('dewp','f4'), ('wndU', 'f4'), ('wndV', 'f4')])
        self.data['pres'] = pres
        self.data['temp'] = temp
        self.data['dewp'] = dewp
        self.data['wndU'] = wndU
        self.data['wndV'] = wndV


        #END NEW INPUT - BEGIN PROGRAM
        #Create figure instance
        fig = plt.figure(figsize=figSize)

        #Create ax2 and make it inverse logarithmic
        self.ax2 = fig.add_subplot(111)
        self.ax2.set_yscale('log')
        self.ax2.set_xscale('linear')
        self.ax2.set_ylim(self.ax2.get_ylim()[::-1])

        #Create ax1 and make it inverse logarithmic and skewed
        self.ax1 = self.ax2.twinx()
        self.ax1.set_yscale('log')
        self.ax1.set_xscale('linear')
        self.ax1.set_ylim(self.ax1.get_ylim()[::-1])
        self.ax1.transLimits = self.ax1.transLimits + Affine2D(np.array([[1.,1.,0.],[0.,1.,0.],[0.,0.,1.]]))
        self.ax1.transData = self.ax1.transScale + (self.ax1.transLimits + self.ax1.transAxes)
        print "Init: ", time.time() - timeCheckStart

        self.new_analyzeParcel(temp, dewp, pres, wndU, wndV)
        
    #Plotting Functions
    def _plotBackground(self):
        timeCheckStart = time.time()
        #Create isobars
        presVals = np.linspace(100,1000,10)
        presVals = np.repeat(presVals,3)
        tempVals = np.array([self.tMin, self.tMax, np.nan])
        tempVals = np.tile(tempVals,(10))
        self.ax2.plot(tempVals, presVals, color=self.colors[2], linestyle=self.lStyle[2], linewidth=self.lWidth[2])
        print "Isob: ", time.time() - timeCheckStart

        timeCheckStart = time.time()
        #Create isotherms
        #(Doing sep calls faster than using np.where() in this case)
        presVals = [self.pMax, self.pMin, np.nan]
        self.ax1.plot([0,0,0], presVals, color=self.colors[4], linestyle=self.lStyle[4], linewidth=self.lWidth[4])
        tempVals = np.linspace(-110,-10,11)
        tempVals = np.repeat(tempVals,3)
        presVals_nE0 = np.tile(presVals,(11))
        self.ax1.plot(tempVals, presVals_nE0, color=self.colors[3], linestyle=self.lStyle[3], linewidth=self.lWidth[3])
        tempVals = np.linspace(10,40,4)
        tempVals = np.repeat(tempVals,3)
        presVals_nE0 = np.tile(presVals,(4))
        self.ax1.plot(tempVals, presVals_nE0, color=self.colors[5], linestyle=self.lStyle[5], linewidth=self.lWidth[5])
        print "Isot: ", time.time() - timeCheckStart

        timeCheckStart = time.time()
        presVals = np.arange(self.pMin, self.pMax+1., 10.)
        #Create dry adiabats
        potTempVals = np.arange(self.tMin + 10., (self.tMax - self.tMin)*2. + self.tMax + 1., 20.) ###Need to come back and make the default maximum automatic.  Would also be nice to have the where statement here again.
        dPresVals = np.append(presVals, np.nan)
        dPresVals = np.tile(dPresVals, (potTempVals.shape[0]))
        potTempVals = np.repeat(potTempVals, presVals.shape[0] + 1)
        tempVals = self._new_getDryAdiabatTemp(potTempVals, dPresVals)
        self.ax1.plot(tempVals, dPresVals, color=self.colors[6], linestyle=self.lStyle[6], linewidth=self.lWidth[6])
        print "DAdi: ", time.time() - timeCheckStart

        timeCheckStart = time.time()
        #Moist Adiabats

        tempVals = np.arange(self.tMin, (self.tMax - self.tMin)*2. + self.tMax + 1., 10.)
        presVals = np.arange(self.pMin, self.pMax+1., 50.) #Note: smaller step values cause issues at lower moist adiabat values near 1000 mb
        #XX NEED TO FIX ABOVE!
        mPresVals = np.append(presVals, np.nan)
        mPresVals = np.tile(mPresVals, (tempVals.shape[0]))
        tempVals = np.repeat(tempVals, presVals.shape[0] + 1) # + 1 to account for nan
        tempVals = self.new_wetLift(1000., tempVals, mPresVals)
        self.ax1.plot(tempVals, mPresVals, color=self.colors[7], linestyle=self.lStyle[7], linewidth=self.lWidth[7])
        
        print "MAdi: ", time.time() - timeCheckStart

        timeCheckStart = time.time()
        #Create mixing ratios
        self.x2TickPres = np.array([self.pMax - (self.pMax - self.mixRatStop)/2.])
        self.w = np.array([1., 2., 3., 5., 8., 13., 21.]) #Note that self.x2TickTemp is dependent upon the number of entries in self.w
        presVals = np.where(presVals > self.mixRatStop, presVals, np.nan)
        self.x2TickTemp = self._new_getMixingRatioTemp(self.w, np.repeat(self.x2TickPres, 7)) 
        for mixRat in self.w:
            tempVals = self._new_getMixingRatioTemp(mixRat, presVals)
            self.ax1.plot(tempVals, presVals, color=self.colors[8], linestyle=self.lStyle[8], linewidth=self.lWidth[8])
        print "MixR: ", time.time() - timeCheckStart

        timeCheckStart = time.time()
        #Create wind profile line
        self.ax2.plot([self.windOffsetPos, self.windOffsetPos], [self.pMax, self.pMin], color=self.colors[9], linestyle=self.lStyle[9], linewidth=self.lWidth[9])
        print "WndP: ", time.time() - timeCheckStart

    def _plotProfile(self):
        timeCheckStart = time.time()
        #Plot temperature, dewpoint, and wind profiles
        self.ax1.plot(self.data['temp'], self.data['pres'], color=self.colors[0], linestyle=self.lStyle[0], linewidth=self.lWidth[0])
        self.ax1.plot(self.data['dewp'], self.data['pres'], color=self.colors[1], linestyle=self.lStyle[1], linewidth=self.lWidth[1])
        self.ax2.barbs(self.windOffsetPos*np.ones(len(self.data['wndU'])), self.data['pres'], self.data['wndU'], self.data['wndV'], color=self.colors[10], linestyle=self.lStyle[10], linewidth=self.lWidth[10])
        print "Plot Profile: ", time.time() - timeCheckStart

    def _plotAxes(self):
        timeCheckStart = time.time()
        #Create temperature labels
        xTickPos = np.arange(self.tMin, self.tMax + 1, 10)
        xTickStr = xTickPos.astype('<U10')
        for t in range(len(xTickStr)):
            xTickStr[t] += u'\u00b0C'
        self.ax1.set_xticks(xTickPos)
        self.ax1.set_xticklabels(xTickStr, clip_on=False)
        print "AxTemp: ", time.time() - timeCheckStart
        
        timeCheckStart = time.time()
        #Create pressure labels ###only plots labels from 1000 - 100 mb! 
        yTickPos = np.arange(100,1001,100)
        yTickStr = yTickPos.astype('<U10')
        for p in range(len(yTickStr)):
            yTickStr[p] += u'mb'
        self.ax2.set_yticks(yTickPos)
        self.ax2.set_yticklabels(yTickStr)
        print "AxPres: ", time.time() - timeCheckStart
        
        timeCheckStart = time.time()
        #Create mixing ratio labels
        for i in range(len(self.x2TickTemp)):
            self.ax1.text(self.x2TickTemp[i], self.x2TickPres[0], str(self.w[i]), ha='center', va='center', alpha=.5, size='smaller')
        print "AxMixR: ", time.time() - timeCheckStart

        #Clears RHS y-axis labels
        self.ax1.set_yticklabels([])

    def _finalizePlot(self):
        #Set plot bounds
        self.ax1.set_ybound(self.pMax, self.pMin)
        self.ax1.set_xbound(self.tMin, self.tMax)
        self.ax2.set_ybound(self.pMax, self.pMin)
        self.ax2.set_xbound(self.tMin, self.tMax)
        
        plt.title(self.title)
        
        #Save and display image
        plt.savefig('skewT.png')
        #plt.show()

    def plot(self):
        #Primary driver function for skewTLogP object
        self._plotProfile()
        self._plotBackground()
        self._plotAxes()
        self._finalizePlot()


    #Helper Functions
    def _new_getMoistAdiabaticLapseRate(self, temp, pres):
        ###NOTE: Must use K and Pa units for inputs.  Don't know why, but it screws up if you don't and convert later
        ###Calculate moist adiabatic lapse rate (see Bluestein)
        sMixR = (1e-3)*self._new_getMixingRatio(self._new_getVaporPres(temp - 273.15), (1/100.)*pres)
        return (278.*temp + sMixR*2.526) / (pres * (1005. + .62198 * sMixR * 2.5e12 / (461.5 * temp**2)))

    def _new_getDryAdiabatTemp(self, potT, pres):
        """
        Calculate temperature along dry adiabat at a given pressure.
        Input:
           potT [C]   : potential temperature
           pres [mb]  : pressure
        Output:
           temp [C]   : temperature along dry adiabat
        Calculation:
            Poisson solution
            temp = ((potT - 273.15) * (pres / 1000)^(gamma)) - 273.15
            where
                gamma = .2858565737 [unitless]
        """
        return (potT + 273.15)*np.power(pres/1000., .2858565737) - 273.15

    def _new_getPotTemp(self, temp, pres):
        """
        Calculate potential temperature for a given temperature and pressure
        Input:
            temp [C]   : temperature
            pres [mb]  : pressure
        Output:
            potT [C]   : potential temperature
        Calculation:
            Poisson solution
            potT = ((temp - 273.15) * (1000 / pres)^(gamma)) - 273.15
            where
                gamma = .2858565737 [unitless]
        """
        return (temp + 273.15)*np.power(1000./pres, .2858565737) - 273.15

    def new_getSatAdiabatTemp(self, sPotT, pres):
        """
        Calculate temperature along a moist adiabat at a given pressure
        Input:
            sPotT [C]   : saturation potential temperature
            pres  [mb]  : pressure
        Output:
            temp  [C]   : temperature
        Calculation:
            (From Stipanuk 1973)
        """
        print ""
        temp = 253.15
        for i in range(12):
            i += 1
            dTemp = (120. / 2**i) * np.sign((sPotT + 273.15) * np.exp(-2.6518986 * self._new_getMixingRatio(self._new_getVaporPres(temp - 273.15), pres) / temp) - (self._new_getPotTemp(temp - 273.15, pres) + 273.15))
            print "A: ", -2.6518986 * self._new_getMixingRatio(self._new_getVaporPres(temp - 273.15), pres) / temp
            print "B: ", self._new_getPotTemp(temp - 273.15, pres) + 273.15
            temp += dTemp
            print "d: ", dTemp
            print "T: ", temp
        print temp - 273.15
        return temp - 273.15
        

    def _new_getSatPotTemp(self, temp, pres):
        """
        Calculate saturation potential temperature for a given temperature and pressure
        Input:
            temp  [C]   : temperature
            pres  [mb]  : pressure
        Output:
            sPotT [C]   : saturation potential temperature
        Calculation:
            (From Stipanuk 1973)
            sPotT = _new_getPotTemp(temp, pres) / exp(b * _new_getMixingRatio(_new_getVaporPres(temp), pres) / (temp + 273.15)) - 273.15
            where
                b = -2.6518986
        """
        return self._new_getPotTemp(temp, pres) / np.exp(-2.6518986 * self._new_getMixingRatio(self._new_getVaporPres(temp), pres) / (temp + 273.15)) - 273.15

    def _new_getMixingRatioTemp(self, mixR, pres):
        """
        Calculate dewpoint (regular) temperature for a given (saturation) mixing ratio and pressure
        Input:
            mixR [g/kg] : (saturation) mixing ratio
            pres [mb]   : pressure
        Output:
            temp [C]  : dewpoint (regular) temperature
        Calculation:
            A = pres * mixR / (6.11 * (622 - mixR))
            B = 1 / 273.15 - 1.846e-4 * ln(A)
            temp = 1 / B - 273.15
        """
        return 1/(1/273.15 - 1.846e-4*np.log(pres*mixR / (6.11*(622 - mixR)))) - 273.15
        #A = np.log10(mixR * pres / (622 + mixR))
        #B = 10 ** (0.0498646455 * A + 2.4082965)
        #C = (10 ** (0.0915 * A) - 1.2035) ** 2
        #return B - 280.23475 + 38.9114 * C

    def _new_getMixingRatio(self, vPres, pres): #W()
        """
        Calculate (saturation) mixing ratio for a given (saturation) vapor pressure and pressure
        Input:
            vPres [mb]   : (saturation) vapor pressure
            pres  [mb]   : pressure
        Output:
            mixR  [g/kg] : (saturation) mixing ratio
        Calculation:
            mixR = 622 * vPres / (100*pres - vPres)
        """
        return 622. * vPres / (pres - vPres) 


    def _new_getVaporPres(self, temp): #ESAT()
        """
        Calculate (saturation) vapor pressure
        Input:
            temp  [C]   : dewpoint (regular) temperature
        Output:
            vPres [mb]  : (saturation) vapor pressure
        Calculation:
            (From Stipanuk 1973 -> Nordquist 1973)
        """
        return 6.11 * np.exp(5417.118093 * (1 / 273.15 - 1 / (temp + 273.15)))
        #temp += 273.15
        #A = 23.832241 - 5.02808 * np.log10(temp)
        #B = 1.3816e7 * 10 ** (11.344 - 0.0303998 * temp)
        #C = 8.1328e3 * 10 ** (3.49149 - 1302.8844 / temp)
        #return 10 ** (A - B + C - 2949.076 / temp)

    def _new_getRH(self, temp, dewp):
        """
        Calculate relative humidity of a parcel
        Input:
            temp [C] : temperature
            dewp [C] : dewpoint temperature
        Output:
            relH [%] : relative humidity
        Calculation:
            relH = 100% * _new_getVaporPres(dewp) / _new_getVaporPres(temp)
        """
        return 100*(self._new_getVaporPres(dewp)/self._new_getVaporPres(temp))

    def _new_getIntersect(tA1, tA2, tB1, tB2, p1, p2):
        #Based on http://paulbourke.net/geometry/lineline2d/
        #HAVEN'T TESTED THIS YET!
        #Assumes pressure levels (y vals) are the same
        uA = ((tB2 - tB1)*(p1 - p2) - (p1 - p2)*(tA1 - tB1)) / ((p1 - p2)*(tA2 - tA1) - (tB2 - tB1)*(p2 - p1))
        return (tA1 + uA*(tA2 - tA1), p1 + uA*(p2 - p1))

    def new_analyzeParcel(self, temp, dewp, pres, wndU, wndV):
        #NOTE: Assumes pressure is decreasing with increasing index
        mixR = self._new_getMixingRatio(self._new_getVaporPres(dewp[0]), pres[0])

        #Calculate LCL (Stipanuk)
        dAdi = self._new_getPotTemp(temp[0], pres[0])
        pLCL = pres[0]
        for i in range(10):
            mRT = self._new_getMixingRatioTemp(mixR, pLCL)
            dAT = self._new_getDryAdiabatTemp(dAdi, pLCL)
            check = 0.02 * (mRT - dAT)
            if abs(check) < .001:
                continue
            pLCL = pLCL * 2 ** check
        tLCL = mRT
        print "---> LCL P, T: ", pLCL, tLCL

        #Calculate Wet Bulb (Builds on LFL calc)
        #sAdi = self._new_getSatPotTemp(tLCL, pLCL)
        #tWbl = self._new_getSatAdiabatTemp(sAdi, pres[0])
        #print "Wet bulb T: ", tWbl

        #Calculate LFC (Builds on Wet Bulb calc)
        #This is an rough calculation. Will likely implement better method / improve this one later.
        #Don't trust this is delta P of measured pressure levels is large!
        #Assumes not saturated at surface
        #sAT = self._new_getSatAdiabatTemp(sAdi, pres)
        #init = np.where(sAT > temp)[0]
        #pLFC = (pres[init[0] - 1] + pres[init[0]]) / 2.
        #tLFC = (temp[init[0] - 1] + temp[init[0]]) / 2.
        #print "---> LFC P, T: ", pLFC, tLFC

        #Calculate EQL (builds on LFC calc)
        #This is an rough calculation. Will likely implement better method / improve this one later.
        #Don't trust this is delta P of measured pressure levels is large!
        #Assumes not saturated at surface
        #A = np.where(pres > pLFC)[0]
        #newP = pres[A[0]:A[-1]]
        #newT = temp[A[0]:A[-1]]
        #sAT = self._new_getSatAdiabatTemp(sAdi, newP)
        #init = np.where(sAT < newT)[0]
        #pEQL = (newP[init[0] - 1] + newP[init[0]]) / 2.
        #tEQL = (newT[init[0] - 1] + newT[init[0]]) / 2.
        #print "---> EQL P, T: ", pEQL, tEQL

        #Calculate CCL
        #This is an rough calculation. Will likely implement better method / improve this one later.
        #Don't trust this is delta P of measured pressure levels is large!
        #Assumes not saturated at surface
        mRT = self._new_getMixingRatioTemp(mixR, pres)
        init = np.where(mRT > temp)[0]
        pCCL = (pres[init[0] - 1] + pres[init[0]]) / 2.
        print "---> CCL P   : ", pCCL
        tCCL = (temp[init[0] - 1] + temp[init[0]]) / 2.

        #Calculate Convective Temp
        dAdi = self._new_getPotTemp(tCCL, pCCL)
        tCon = self._new_getDryAdiabatTemp(dAdi, pres[0])
        print "---> Conv T  : ", tCon

        #Plot info on skew t
        self.ax1.text(tLCL, pLCL, "  LCL", ha='left', va='center', size='smaller')
        self.ax1.plot(tLCL, pLCL, 'k_', markeredgewidth=2, markersize=11)
        self.ax1.text(tCCL, pCCL, "  CCL", ha='left', va='center', size='smaller')
        self.ax1.plot(tCCL, pCCL, 'k_', markeredgewidth=2, markersize=11)
        self.ax1.text(tCon, pres[0], " TCon", ha='left', va='center', size='smaller')
        self.ax1.plot(tCon, pres[0], 'kx', markeredgewidth=2)

#####
#####
    #This bit based on functions found in nsharp's thermo.c
    def wobusFunc(self, temp):
    #/*  pres             - Pressure to raise parcel (mb)         */
    #/*  thm              - Sat. Pot. Temperature of parcel (c)   */

        temp = temp - 20
        if (temp <= 0):
            A = 1. + temp*(-8.841660499999999e-03 + temp*(1.4714143e-04 + temp*(-9.671989000000001e-07 + temp*(-3.2607217e-08 + temp*(-3.8598073e-10)))))
            return 15.13 / A**4
        else:
            A = temp*(4.9618922e-07 + temp*(-6.1059365e-09 + temp*(3.9401551e-11 + temp*(-1.2588129e-13 + temp*(1.6688280e-16)))))
            A = 1 + temp*(3.6182989e-03 + temp*(-1.3603273e-05 + A))
            return 29.94 / A**4 + .96*temp - 14.8

    def new_wobusFunc(self, temp):
        #Wobus function designed to work with an array of temps
        #Used in conjunction with new_wobusFunc_lt0 and new_wobusFunc_gt0
        temp = temp - 20
        return np.where(temp <= 0, self.new_wobusFunc_lt0(temp), self.new_wobusFunc_gt0(temp))
    def new_wobusFunc_lt0(self, temp):
            A = 1. + temp*(-8.841660499999999e-03 + temp*(1.4714143e-04 + temp*(-9.671989000000001e-07 + temp*(-3.2607217e-08 + temp*(-3.8598073e-10)))))
            return 15.13 / A**4
    def new_wobusFunc_gt0(self, temp):
            A = temp*(4.9618922e-07 + temp*(-6.1059365e-09 + temp*(3.9401551e-11 + temp*(-1.2588129e-13 + temp*(1.6688280e-16)))))
            A = 1 + temp*(3.6182989e-03 + temp*(-1.3603273e-05 + A))
            return 29.94 / A**4 + .96*temp - 14.8

    def wetLift(self, presStart, tempStart, presEnd):
        #p       (float)         Pressure of initial parcel (hPa)
        #t       (float)         Temperature of initial parcel (C)
        #p2      (float)         Pressure of final level (hPa)

        # =>Temperature (C [float])
        potT = self._new_getPotTemp(tempStart, presStart)
        corrT = potT - self.wobusFunc(potT) + self.wobusFunc(tempStart)
        return self.satLift(presEnd, corrT)
    def new_wetLift(self, presStart, tempStart, presEnd):
        #p       (float)         Pressure of initial parcel (hPa)
        #t       (float)         Temperature of initial parcel (C)
        #p2      (float)         Pressure of final level (hPa)

        # =>Temperature (C [float])
        potT = self._new_getPotTemp(tempStart, presStart)
        corrT = potT - self.new_wobusFunc(potT) + self.new_wobusFunc(tempStart)
        return self.new_satLift(presEnd, corrT)

    def satLift(self, pres, temp):
        #/*  Returns the temperature (c) of a parcel (thm),           */
        #/*  when lifted to level (pres).                             */
        #/*                                                           */
        #/*  pres             - Pressure to raise parcel (mb)         */
        #/*  thm              - Sat. Pot. Temperature of parcel (c)   */

        if np.abs(pres - 1000.) - 1.e-3 <= 0:
            return temp

        e0 = 999.
        while (np.abs(e0) - .1 > 0):
            if e0 == 999.:
                powF = np.power(pres/1000., .2858565737)
                t1 = self._new_getPotTemp(temp, pres)
                woto = self.wobusFunc(t1)
                wotm = self.wobusFunc(temp)
                e1 = woto - wotm
                rate = 1.
            else:
                rate = (t2 - t1) / (e2 - e1)
                t1 = t2
                e1 = e2
            t2 = t1 - e1 * rate
            e2 = (t2 + 273.15) / powF - 273.15
            wot2 = self.wobusFunc(t2)
            woe2 = self.wobusFunc(e2)
            e2 = e2 + wot2 - woe2 - temp
            e0 = e2 * rate

        return t2 - e0

    def new_satLift(self,pres,temp):
        return np.where(np.abs(pres - 1000.) <= 1.e-3, temp, self.new_satLift_gtVal(pres, temp))
    def new_satLift_gtVal(self,pres,temp):
        e0 = 999.* np.ones(temp.shape)
        start = 'Y'
        #left = temp.shape[0]
        count = 1
        while np.any(np.abs(e0) - .1 > 0):
            if count >= 100:
                print("ERROR: Call to function new_satLift_gtVal() failed after {0} iterations.").format(count)
            notDoneIdx = np.where(np.abs(e0) > .1)[0]
            count += 1
            if start == 'Y':
                start = 'N'
                powF = np.power(pres/1000., .2858565737)
                t1 = self._new_getPotTemp(temp, pres)
                woto = self.new_wobusFunc(t1)
                wotm = self.new_wobusFunc(temp)
                e1 = woto - wotm
                rate = np.ones(temp.shape)
                #Just initialize these values as arrays - they will be overwritten
                t2 = np.ones(temp.shape)
                e2 = np.ones(temp.shape)
                wot2 = np.ones(temp.shape)
                woe2 = np.ones(temp.shape)
                
            else:
                rate[notDoneIdx] = (t2[notDoneIdx] - t1[notDoneIdx]) / (e2[notDoneIdx] - e1[notDoneIdx])
                t1[notDoneIdx] = t2[notDoneIdx]
                e1[notDoneIdx] = e2[notDoneIdx]
            #Values correct until here
            t2[notDoneIdx] = t1[notDoneIdx] - e1[notDoneIdx] * rate[notDoneIdx]
            e2[notDoneIdx] = (t2[notDoneIdx] + 273.15) / powF[notDoneIdx] - 273.15
            wot2[notDoneIdx] = self.new_wobusFunc(t2[notDoneIdx])
            woe2[notDoneIdx] = self.new_wobusFunc(e2[notDoneIdx])
            e2[notDoneIdx] = e2[notDoneIdx] + wot2[notDoneIdx] - woe2[notDoneIdx] - temp[notDoneIdx]
            e0[notDoneIdx] = e2[notDoneIdx] * rate[notDoneIdx]

        return t2 - e0

#####

#####

x = skewTLogP()
x.plot()
#x.solve(0., 1000, 1.2e-3)
#x.new_getSatAdiabatTemp(np.array([-40, -20,-10,0,10,20,40]), np.array([1000,1000,1000,1000,1000,1000,1000]))
#p = np.array([1000.,900.,800.,700.,600.,500.,400.,300.,200.,100.])
#t = np.array([20.  ,20., 20., 20., 20., 20., 20., 20., 20., 20.])
#tk = 273.15 + np.array([20.  ,20., 20., 20., 20., 20., 20., 20., 20., 20.])
#print x.TSA(np.array([273.15]), np.array([1000.]))
#print x.OS(np.array([273.15]), np.array([1000.]))

endTime = time.time()
totalTime = endTime - startTime
print("Total time: {0} s".format(totalTime))
