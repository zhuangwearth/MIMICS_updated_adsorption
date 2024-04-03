#' MIMICSb: Implementation of the MIMICS model with revised mineral adsorption description.
#' @references Wieder et al. 2015 GCB.
#' The codes of original version are available at https://github.com/wwieder/MIMICS/releases/tag/MIMICS_v0.1 
#' The following version was revised by Wei Zhuang in Dec 2022
#' Author's email: zhuangwearth@163.com

#' Main improvement
# (1) Compared with MIMICSb,  the maximum adsorption capacity is revised hereas:
#     Qmax = 86*(ClAY + SILT)/100

#' Note:
# the codes are for bare follow experiments. Thus, litter pools are excluded. 
# There are five pools:
# MIC_1 is r-type microbial biomass
# MIC_2 is r-type microbial biomass
# SOM_1 is physically protected SOM, i.e. PHYS SOM for short
# SOM_2 is chemically recalcitrant SOM, i.e. CHEM SOM for short
# SOM_3 is available SOM, i.e. AVAIL SOM for short

MIMICS_BF <- function 
( beta,               # optimized, microbial biomass density dependence factor
  t,                  # date index
  timestep,           # equal to 1, for hourly step; 24, daily step; 24*365, yearly step
  TSOI,               # soil temperature (degree)
  CLAY,               # Clay content (%)
  SILT,               # Silt content (%)
  LIG,                # ligin content  
  CN,                 # the ratio of C to N
  bulkD,              # bulk density
  
  # the values of the following parameter except for Tao_MOD and r_ads_des
  # are adopted from Wieder et al. (2015)
  Vslope    = 0.063,
  Vint      = 5.47,
  aV        = 8e-6,
  Kslope_SOMa = 0.017,
  Kint      = 3.19,
  aK        = 10,
  Tao_MOD   = 1.0,
  CUE_r_SOMa= 0.55,
  CUE_k_SOMa= 0.75,

  fCLAY     = CLAY/100,                    # clay content, convert from clay fraction to %
  calCN     = (1 / CN) / 2.5 * 100, 
  fMET      = 0.85 - 0.013 * LIG /calCN,             #as calculated in DAYCENT
  
  # Calculate Vmax(T) & (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
  Vslopes   = array(Vslope,dim=6),
  Vmax      = exp(TSOI * Vslopes + Vint) * aV,
  MOD1      = c(10, 2, 10, 3, 3, 2), 
  VMAX      = Vmax * MOD1  
              * timestep,    
  
  #Calculate Km(T)    
  # Km = f(x = TSOI, pars = c(Kslope, Kint, aK, MOD2, pSCALAR)) 
  Kslopes    = c(0.017, 0.027, Kslope_SOMa, 0.017, 0.027, Kslope_SOMa),
  Km        = exp(Kslopes * TSOI + Kint) * aK,
  
  pSCALAR   = 1.0/ (2.0 * exp(-2.0*(sqrt(fCLAY)))), #Scalar for texture effects on PHYS SOM
  MOD2      = c(0.125, 0.5,  0.25 * pSCALAR, 
                0.5,   0.25, pSCALAR/6),   
  KM        = Km * MOD2,
  
  # parameter related to mic turnover --------------
  Tao_MOD   = 1.0,   # Tao_MOD is set as 1.0 in the bare fallow experiments.   
  tao       = c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))* Tao_MOD  # microbial biomass turover rates of r and K types 
                * timestep,       
  
  # The microbial residues will be partitioned to:
  fPHYS     = 0,     # fraction to PHYS SOM is set zero in the version.
  fCHEM     = c(0.1 * exp(-3*fMET)  , 0.3 * exp(-3*fMET)  ),   # fraction to CHEM SOM
  fAVAI     = 1- (fPHYS + fCHEM),
  
  # desorb
  desorb    = 1.5e-5 * exp(-1.5*(fCLAY))      #CHANGED FOR GLOBAL RUN!!!   
              *timestep,
  
  # Adsorption 
  r_ads_des = 6,                  # ratio of adsorption rate to desorption flux
  adsorb    = r_ads_des * desorb, # Huang et al (2018) 
  
  # Qmax, the maximum adsorption capacity
  Qmax_v    = 86*(CLAY + SILT)/100,         # mg C per gSoil, most soils are in HM class
  Qmax      = Qmax_v * bulkD,                              # mg C per cm3 Soil

  # cue
  CUE       = c(CUE_r_SOMa, 0.25, CUE_k_SOMa, 0.35),
  KO        = c(4, 4), 
  
  ival=c(MIC_1=0.1, MIC_2=0.1, SOM_1=8, SOM_2=5, SOM_3=3) # the initial state of C pools
)
{
  t_start=min(t)
  t_end=max(t)
  nr=5	
  f=function(C,t){
    MIC_1 = C[1]
    MIC_2 = C[2] 
    SOM_1 = C[3] # physical SOM
    SOM_2 = C[4] # chemical SOM
    SOM_3 = C[5]	  
    O=matrix(byrow=TRUE,nrow=nr,c(MIC_1^beta * tao[1], # Modified
                                  MIC_2^beta * tao[2], # Modified
                                  desorb * SOM_1,      # as original MIMICS
                                  ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) + 
                                   (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2))),  # oxidation of chemically recalcitrant SOM 
                                  (MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3) +    # decomp of AVAIL SOM by MIC_1
                                   MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3) +    # decomp of AVAIL SOM by MIC_2
                                   SOM_3 * adsorb * (1 - SOM_1/Qmax) )
    )
    )
    return(O)
  }
  alpha=list()
  ##Percentage of Flows from MIC_1 (r type) 
  alpha[["1_to_3"]]=function(C,t){
    fPHYS[1] # to PHYS SOM
  }
  alpha[["1_to_4"]]=function(C,t){
    fCHEM[1] # to CHEM SOM
  }
  alpha[["1_to_5"]]=function(C,t){
    fAVAI[1] # to AVAI SOM
  }
  
  ##Percentage of Flows from MIC_2 (K type) 	
  alpha[["2_to_3"]]=function(C,t){
    fPHYS[2] # to PHYS SOM
  }
  alpha[["2_to_4"]]=function(C,t){
    fCHEM[2] # to CHEM SOM
  }
  alpha[["2_to_5"]]=function(C,t){
    fAVAI[2] # to AVAI SOM
  }
  
  ##Percentage of Flows from SOM_1 (physical SOM) 
  alpha[["3_to_5"]]=function(C,t){
    1 # all the lost AVAI SOM go to here
  }
  
  ##Percentage of Flows from SOM_2 (CHEM SOM) 
  alpha[["4_to_5"]]=function(C,t){
    1 # all the lost CHEM SOM go to here
  }	
  
  ##Percentage of Flows from SOM_3 (Available SOM) 
  # to MIC_1  
  alpha[["5_to_1"]]=function(C,t){
    CUE[1] * C[1] * VMAX[3] * C[5] / (KM[3] + C[5]) /
	(C[1] * VMAX[3] * C[5] / (KM[3] + C[5]) + C[2] * VMAX[6] * C[5] / (KM[6] + C[5]) + 
	 C[5] * adsorb * (1 - C[3]/Qmax) )
  } 
  # to MIC_2 
  alpha[["5_to_2"]]=function(C,t){  
    CUE[3] * C[2] * VMAX[6] * C[5] / (KM[6] + C[5]) /
	(C[1] * VMAX[3] * C[5] / (KM[3] + C[5]) + C[2] * VMAX[6] * C[5] / (KM[6] + C[5]) +
	 C[5] * adsorb * (1 - C[3]/Qmax) )
  }	
  # to PHYS SOM (absorption)
  alpha[["5_to_3"]]=function(C,t){  
    C[5] * adsorb * (1 - C[3]/Qmax)  /
	(C[1] * VMAX[3] * C[5] / (KM[3] + C[5]) + C[2] * VMAX[6] * C[5] / (KM[6] + C[5]) +
	 C[5] * adsorb * (1 - C[3]/Qmax) )
  }	
  
  Anl=new("TransportDecompositionOperator",t_start,Inf,nr,alpha,f)
  inputrates=BoundInFluxes(
    function(t){
      matrix(
        nrow=nr,
        ncol=1,
        c(0,0, 0, 0, 0)
      )
    },
    t_start,
    t_end
  )
  modnl=GeneralNlModel(t, Anl, ival, inputrates, deSolve.lsoda.wrapper)
  return(modnl)
}
