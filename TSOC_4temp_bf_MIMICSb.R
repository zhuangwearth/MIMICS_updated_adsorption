# model to get total carbon (TSOC)
# # the initial total carbon has an unit of mgC_per_cm3

# calculate TSOC for both bare fallow and incubation
TSOC_4temp_bf_MIMICS = function(driver, pars){
	    
        fMIC_1 = pars[1]; fMIC_2 = pars[2]
        fSOM_1 = pars[3]; fSOM_3 = pars[4]
        fSOM_2 = 1- fMIC_1 - fMIC_2 - fSOM_1 - fSOM_3
        beta   = pars[5]
        
		# for bare fallow
		# Site	Askov	Grignon	   Ultuna	Versailles
        # MAT	7.8	     10.7		5.5		10.7
        
		ta = MAT
		exp_type = 'LTBF'
		df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
		cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
		
		Ex= MIMICS_BF(t = df_driver$date_sin,  # date since the start of bare fallow, years
                      timestep = 24*365, # yearly
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= fMIC_1 * cint, MIC_2=fMIC_2*cint, 
                               SOM_1= fSOM_1 * cint, SOM_2=fSOM_2*cint, SOM_3= fSOM_3*cint))
        Ct=getC(Ex)
        TSOC_pre_bf =rowSums(Ct) # total carbon storage	
		
		# Get final size of carbon pools
        pools = as.data.frame(cbind(df_driver$date_sin, Ct, TSOC_pre_bf))
		colnames(pools) = c('date_sin', 'MIC_1', 'MIC_2', 'SOM_1', 'SOM_2', 'SOM_3','TSOC')
        TSOC_fin = pools[pools$date_sin == last,]
        MIC_1fin = TSOC_fin$MIC_1
        MIC_2fin = TSOC_fin$MIC_2
        SOM_1fin = TSOC_fin$SOM_1
        SOM_2fin = TSOC_fin$SOM_2
        SOM_3fin = TSOC_fin$SOM_3
		
		
		## for incubation
		## Temperature 4 for initial sample  
        ta = 4 
        exp_type = 'initial.inc'
        df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
        cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
        Ex= MIMICS_BF(t = df_driver$date_sin, 
                      timestep = 24, # daily
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= fMIC_1 * cint, MIC_2=fMIC_2*cint, 
                               SOM_1= fSOM_1 * cint, SOM_2=fSOM_2*cint, SOM_3= fSOM_3*cint))
        Ct=getC(Ex)
        TSOC_pre_4ini =rowSums(Ct) # total carbon storage

		## Temperature 4 for final sample  
        exp_type = 'final.inc'
        df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
        cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
        Ex= MIMICS_BF(t = df_driver$date_sin, 
                      timestep = 24, # daily
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= MIC_1fin, MIC_2= MIC_2fin, 
                               SOM_1= SOM_1fin, SOM_2= SOM_2fin, SOM_3= SOM_3fin))
        Ct=getC(Ex)
        TSOC_pre_4fin =rowSums(Ct) # total carbon storage


		## Temperature 12 for initial sample  
        ta = 12
        exp_type = 'initial.inc'
        df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
        cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
        Ex= MIMICS_BF(t = df_driver$date_sin, 
                      timestep = 24, # daily
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= fMIC_1 * cint, MIC_2=fMIC_2*cint, 
                               SOM_1= fSOM_1 * cint, SOM_2=fSOM_2*cint, SOM_3= fSOM_3*cint))
        Ct=getC(Ex)
        TSOC_pre_12ini =rowSums(Ct) # total carbon storage

		## Temperature 12 for final sample  
        exp_type = 'final.inc'
        df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
        cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
        Ex= MIMICS_BF(t = df_driver$date_sin, 
                      timestep = 24, # daily
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= MIC_1fin, MIC_2= MIC_2fin, 
                               SOM_1= SOM_1fin, SOM_2= SOM_2fin, SOM_3= SOM_3fin))
        Ct=getC(Ex)
        TSOC_pre_12fin =rowSums(Ct) # total carbon storage


		## Temperature 20 for initial sample  
        ta = 20
        exp_type = 'initial.inc'
        df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
        cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
        Ex= MIMICS_BF(t = df_driver$date_sin, 
                      timestep = 24, # daily
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= fMIC_1 * cint, MIC_2=fMIC_2*cint, 
                               SOM_1= fSOM_1 * cint, SOM_2=fSOM_2*cint, SOM_3= fSOM_3*cint))
        Ct=getC(Ex)
        TSOC_pre_20ini =rowSums(Ct) # total carbon storage

		## Temperature 20 for final sample  
        exp_type = 'final.inc'
        df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
        cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
        Ex= MIMICS_BF(t = df_driver$date_sin, 
                      timestep = 24, # daily
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= MIC_1fin, MIC_2= MIC_2fin, 
                               SOM_1= SOM_1fin, SOM_2= SOM_2fin, SOM_3= SOM_3fin))
        Ct=getC(Ex)
        TSOC_pre_20fin =rowSums(Ct) # total carbon storage


		## Temperature 35 for initial sample  
        ta = 35
        exp_type = 'initial.inc'
        df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
        cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
        Ex= MIMICS_BF(t = df_driver$date_sin, 
                      timestep = 24, # daily
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= fMIC_1 * cint, MIC_2=fMIC_2*cint, 
                               SOM_1= fSOM_1 * cint, SOM_2=fSOM_2*cint, SOM_3= fSOM_3*cint))
        Ct=getC(Ex)
        TSOC_pre_35ini =rowSums(Ct) # total carbon storage

		## Temperature 4 for final sample  
        exp_type = 'final.inc'
        df_driver = driver[driver$temp == ta &  driver$exp_type == exp_type,]
        cint = df_driver[1, 'SOC_mgC_per_cm3']   # the initial total carbon
        Ex= MIMICS_BF(t = df_driver$date_sin, 
                      timestep = 24, # daily
                      bulkD = ibulkD,
                      TSOI = ta,
                      CLAY = CLAY,
                      beta = beta,
                      ival = c(MIC_1= MIC_1fin, MIC_2= MIC_2fin, 
                               SOM_1= SOM_1fin, SOM_2= SOM_2fin, SOM_3= SOM_3fin))
        Ct=getC(Ex)
        TSOC_pre_35fin =rowSums(Ct) # total carbon storage

        TSOC_pre = c(TSOC_pre_bf, TSOC_pre_4ini,  TSOC_pre_12ini,  TSOC_pre_20ini,  TSOC_pre_35ini,
                                  TSOC_pre_4fin,  TSOC_pre_12fin,  TSOC_pre_20fin,  TSOC_pre_35fin)
        }
