#' MIMICSa: Implementation of the MIMICS model with default mineral adsorption description.
#' @references Wieder et al. 2015 GCB.
#' The codes of original version are available at https://github.com/wwieder/MIMICS/releases/tag/MIMICS_v0.1 
#' The following version was revised by Wei Zhuang. in Dec 2022
#' Author's email: zhuangwearth@163.com

#' Main improvement
# (1) The model reconstructed in the framework of SoilR package by Sierra et al. (2012)
# (2) microbial biomass density dependence factor (beta) was added

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
  LIG,                # ligin content  
  CN,                 # the ratio of C to N
  
  # the values of the following parameters except for Tao_MOD
  # are adopted from Wieder et al. (2015)
  Vslope    = 0.063,
  Vint      = 5.47,
  aV        = 8e-6,
  Kslope_SOMa = 0.017,
  Kint      = 3.19,
  aK        = 10,

  CUE_r_SOMa= 0.55,
  CUE_k_SOMa= 0.75,
  
  fCLAY     = CLAY/100,                    # clay content, convert from clay fraction to %
  calCN     = (1 / CN) / 2.5 * 100, 
  fMET      = 0.85 - 0.013 * LIG /calCN,             #as calculated in DAYCENT

  # Calculate Vmax(T) & (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
  Vslopes    = array(Vslope,dim=6),
  Vmax      = exp(TSOI * Vslopes + Vint) * aV,
  MOD1      = c(10, 2, 10, 3, 3, 2), 
  VMAX      = Vmax * MOD1  
              *timestep,    # convert date unit
  
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
  tao       = c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))* Tao_MOD  # microbial biomass turnover rates of r and K types 
                *timestep,       # convert date unit
  
  # The microbial residues will be partitioned to:
  fPHYS     = c(0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY)),   #fraction to PHYS SOM
  fCHEM     = c(0.1 * exp(-3*fMET)  , 0.3 * exp(-3*fMET)  ),   #fraction to CHEM SOM
  fAVAI     = 1- (fPHYS + fCHEM),
  
  # desorption rate
  desorb    = 1.5e-5 * exp(-1.5*(fCLAY))      #CHANGED FOR GLOBAL RUN!!!   
              *timestep,    # convert unit
  # cue
  CUE       = c(CUE_r_SOMa, 0.25, CUE_k_SOMa, 0.35),
  KO        = c(4, 4), 
  
  ival=c(MIC_1=0.1, MIC_2=0.1, SOM_1=8, SOM_2=5, SOM_3=3) # the initial state of C pools
)
{
  t_start=min(t)
  t_end=max(t)
  nr=5	# five pools
  f=function(C,t){
    MIC_1 = C[1]
    MIC_2 = C[2] 
    SOM_1 = C[3] 
    SOM_2 = C[4] 
    SOM_3 = C[5]	  
    O=matrix(byrow=TRUE,nrow=nr,c(MIC_1^beta * tao[1],  # Modified
                                  MIC_2^beta * tao[2],  # Modified
                                  SOM_1 * desorb,       # desorption flux 
                                  ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) + 
                                     (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2))),  # oxidation of chemically recalcitrant SOM 
                                  MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3) +             # decomp of AVAIL SOM by MIC_1
                                    MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)             # decomp of AVAIL SOM by MIC_2
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
  
  ##PerCentage of Flows from SOM_1 (physical SOM) 
  alpha[["3_to_5"]]=function(C,t){
    1 # all the lost AVAI SOM go to here
  }
  
  ##Perentage of Flows from SOM_2 (CHEM SOM) 
  alpha[["4_to_5"]]=function(C,t){
    1 # all the lost CHEM SOM go to here
  }	
  
  ##Perentage of Flows from SOM_3 (Available SOM) 
  # to MIC_1  
  alpha[["5_to_1"]]=function(C,t){
    CUE[1] * C[1] * VMAX[3] * C[5] / (KM[3] + C[5]) /
	(C[1] * VMAX[3] * C[5] / (KM[3] + C[5]) + C[2] * VMAX[6] * C[5] / (KM[6] + C[5]))
  }
  
  ##Perentage of Flows from SOM_3 (Available SOM) 
  # to MIC_2 
  alpha[["5_to_2"]]=function(C,t){  
    CUE[3] * C[2] * VMAX[6] * C[5] / (KM[6] + C[5]) /
	(C[1] * VMAX[3] * C[5] / (KM[3] + C[5]) + C[2] * VMAX[6] * C[5] / (KM[6] + C[5]))
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
