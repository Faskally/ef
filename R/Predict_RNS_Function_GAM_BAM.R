# ----------------------------------------------------------------------------------------
# Function to allow you to predict with and without river network smoothers in gam/bam models
# ----------------------------------------------------------------------------------------

# This function allows you to predict for data in the following formats, using gam or bam models:
  # With river networks smoothers where all prediction data points have an RNS 
    # e.g. predicting for observed catchments or lifestages
  # Without river networks smoothers where no prediction data points have an RNS 
    # e.g. predicting to observed catchments or lifestages
  # With river networks smoothers AND without river network smoothers in the 
  # prediction dataset 
    # e.g. predicting to both observed AND unobserved catchments or lifestages

# The function has been optimised so there is no longer a requirement to subset data,
# predictions for all data can be done in one step. Furthermore, there is no need
# to provide a vector of RNS names or random effects names within the function call 
# these are now taken from the model object and the column used to generate RNS or RE

# The function can also deal with circumstances where the RNS may not have multiple levels 
# e.g. a model for a single catchment or lifestage. Predictions can be made with and 
# without RNS as above

# Required libraries are mgcv and ef (although these are called within the function so 
# providing they are available on your computer they will be loaded by the function)

# ---------------
# Function call specification: 
# ---------------

# predict_RNS_gam_bam <- function(gmod, newdata, RNS_col, by_factor_col,
                                #to_drop_col=NULL, type = "link", se.fit = FALSE) 

# gmod = Model that you are predicting with

# newdata = Prediction dataset (data points you are predicting for)

# RNS_col = column for the river network smoother 'riversmooth' in our case

# by_factor_col = column/factors that the river network smoother is different for 
                  # (e.g. by= in model call) in our case this may be CTMName or Lifestage

# to_drop_col = the default is NULL - here any column names you wish to drop can be 
                # provide e.g. if you wanted to predict without a hydrometric area smoother
                # NOTE: Random effects columns are automatically removed from the lpmatrix 
                # within the function and no longer need to be provided.Similarly for catchments
                # without an RNS you do not need to provide the RNS columns. If you want to predict
                # without RNS, say for partial effects plots, just add in the RNS column
                # (in our case 'riversmooth) and function will generate the appropriate names to drop
                # from lpmatrix

# single_RNS = the default is FALSE - RNS will have multiple levels. If the RNS has only one level
                # e.g. a single catchment or lifestage this should be changed to true

# single_RNS_factor_col = the default is NULL - if the RNS has only one level (e.g. a single 
                # catchment or lifestage) the column where this level could be found in a prediction
                # data should be included here ('CTMName' or 'Lifestage' in our case)

# single_RNS_factor = the default is NULL - if the RNS has only one level (e.g. a single 
                # catchment or lifestage) the level in the model should be included here 
                # (e.g. "River Dee at mouth" or "Fry" in our case)

# type =  the default is "link" (predictions) but can use lpmatrix

# se.fit = the default is FALSE (no standard errors)
          # NOTE: if se.fit=T the result returned is a list but using $fit or $se.fit after
          # the model call lets you add each in a column

# Notes: 03/09/2018 - FLJ


predict_RNS_gam_bam <- function(gmod, newdata, RNS_col, by_factor_col, to_drop_col=NULL, 
                                  single_RNS = FALSE, single_RNS_factor_col = NULL,
                                  single_RNS_factor = NULL,
                                  type = "link", se.fit = FALSE) {
    # load neccessary libraries
    require(mgcv)
    require(ef)
    #browser() # for troubleshooting
    
    # ---------------------------------------------------
    # Generate information required from the model object
    # ---------------------------------------------------
    
    # Get original order of the prediction dataset which you can use to make sure your ordering 
    # stays consistent when data is subsetted and joined back together
    org_order<-rownames(newdata)
    
    # --------
    # Extract random effects from model so don't need to specify 
    # them as a vector in the function call
    # --------
    # Return the position of model terms which are random effects
    re.position<-lapply(gmod$smooth, function(x) attributes(x)$class[1] == "random.effect")
    
    # Extract the smoother infomation for each random effect
    random.effects<-gmod$smooth[grep(TRUE, re.position)]
    
    # Generate the names of all the columns in the lpmatrix that are the random effects,
    # which you want to remove so they are not predicted with 
    # Generate the names by getting all the point values, as depending on the k values 
    # (/ N of data points for random effect) the random effect name will be .1 to .(length(RE levels))
    random.effects.names <- unlist(lapply(random.effects, function(x) paste0(x$label, ".", 1:x$rank)))
    
    # NOW READY TO BEGIN PREDICTIONS BASED ON RULES ESTABLISHED IN THE FUNCTION CALL
    
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # Check if model has an RNS which varies by factor level using single_RNS == TRUE
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
  
    # IF RNS HAS ONLY ONE LEVEL (E.G. SINGLE CATCHMENT OR LIFESTAGE DO THE FOLLOWING:
    if(single_RNS == TRUE){
      cat("Single RNS level e.g. single CTM model")
      
      # ------------------------------------  
      # Subset to data by if you will be predicting with and without the RNS 
      # ------------------------------------  
      
      # Get data with RNS 
      newdata_RNS<-newdata[as.character(newdata[, c(single_RNS_factor_col)]) %in% 
                             single_RNS_factor,]
      
      # Get data without RNS
      # This could end up being empty if you are not make predictions where you don't have RNS 
      # e.g. only predicting in catchments with data or only predicting for lifestages the model 
      # was fitted to
      newdata_noRNS<-newdata[! as.character(newdata[, c(single_RNS_factor_col)]) %in% 
                               single_RNS_factor,]
      
      # Next use a series of if statements which determine if you are predicting only with the RNS, 
      # only without RNS or with AND without RNS as each have to be generated slightly differently
      
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # GET PREDICTIONS WITH & WITHOUT RIVER NETWORK SMOOTHERS 
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      
      if(nrow(newdata_RNS) > 0 & nrow(newdata_noRNS) > 0) {
        cat("Predictions with and without RNS")  
        
        # ---------------------------------------------------------------
        # Generate predictions with RNS
        # ---------------------------------------------------------------
        
        # Generate the lpmatrix 
        # We know that the RNS bits of this will be wrong and thus we need to update them 
        X <- predict(gmod, newdata = newdata_RNS, type = "lpmatrix")
        
        # Find appropriate bits of the smooth parts of model (held in gmod$smooth)
        # that are RNS bits using the relevant column that refers to the RNS -in our case 'riversmooth'
        which.smooth <- grep(RNS_col, sapply(gmod$smooth, "[[", "term")) 
        
        # Extract the smoother info for each RNS 
        # This will be a list if the RNS is by a factor (as each are held in a separate list element)
        smooth<-gmod$smooth[which.smooth]
        
        # Get names of all the RNS coefficients
        # The point values (pulled from the rank) relate to the k value provided e.g. 
        # if RNS has K value of 3 (3 df) the RNS will be RNS:FACTOR.1, RNS:FACTOR.2 RNS:FACTOR.3
        coef.names <- unlist(lapply(smooth, function(x) paste0(x$label, ".", 1:x$rank)))
        
        # Find all the RNS covariate names i.e. bits of data used as RNS (
        # in our case Order_Smooth_1 Order_Smooth_2 etc) and reduce the column names 
        # to the number of columns that can be used in the model i.e. to the K value (pulled from the rank)
        # e.g. if RNS has K value of 3 (3 df) the RNS will be Order_Smooth_1, Order_Smooth_2 Order_Smooth_3
        covar.names <-unlist(lapply(smooth, function(x) colnames(x$xt$X)[1:x$rank + ncol(x$xt$X)- x$rank]))
        
        # Replace the elements in the lpmatrix relating to the RNS bases with the appropriate values 
        X[,coef.names] <- as.matrix(newdata_RNS[, covar.names])
        
        # Keep copy for lmpatrix export which has correct RNS columns but no columns removed
        X_full<-X
        
        # Find which bits of the lpmatrix you want to drop (i.e. not predict with) 
        # Based on the columns provided in the function call and any random effects 
        # (determined from the model at the start) 
        if(is.null(to_drop_col)){
          # If nothing is provide in the to_drop_col just list the random effects
          drop_names<-random.effects.names
        } else {
          # If to_drop_col is provided create a vector of the columns in the lpmatrix which
          # refer to that column (collapse "|" allows a vector of multiple columns to drop
          # to be provided in the function call as these will be looked for separately) and 
          # any random effects
          drop_names<-c(random.effects.names, grep(paste(to_drop_col,collapse="|"), colnames(X), value=T))
        }
        
        # Remove any columns you do not want to predict with from the lpmatrix 
        # and model coefficients
        if(!is.null(drop_names)){
          # Remove all the lpmatrix you don't want to predict with i.e. the random effects
          X <- X[, ! colnames(X) %in% drop_names]
          # Remove all the coefficients you don't want to predict with i.e. the random effects
          coefs <- coef(gmod)[setdiff(names(coef(gmod)), drop_names)]
        } else {
          # If drop_names is null that means the model has no random effects and no columns which 
          # you do not want to predict with, so you do not need to update the lpmatrix or coefficients 
          coefs <- coef(gmod)
        }
        
        # ------------------------------
        # GENERATE THE RESULTS ASKED FOR IN THE FUNCTION CALL 
        # ------------------------------
        if (type == "lpmatrix") {
          cat("RNS lpmatrix generated") # Don't need to do anything lpmatrix is generated already
        } else if (type == "link") {
          # Generate the fit but also give it the names of the original 
          # row order (from the lpmatrix rownames)
          fit <- setNames(as.vector(X %*% coefs), rownames(X))
          
          if (se.fit) {# if se.fit=FALSE (default) just return fit 
            # If se.fit is true get the standard errors too:
            if(!is.null(drop_names)){
            # Remove any columns you do not want to predict with from the lpmatrix 
            # and model coefficients
            # Get positions of the columns/rows we don't want want to predict with
            drop_names_pos<-which(names(gmod$coefficients) %in% drop_names)
            # Remove columns in Vp (covariance matrix) which we don't want to predict with 
            Vp <- gmod$Vp[-drop_names_pos, -drop_names_pos]
            } else {
              # If drop_names is null that means the model has no random effects and no columns which 
              # you do not want to predict with, so you do not need to update the lpmatrix or coefficients 
            Vp <- gmod$Vp 
            }
            # Get standard errors row by row to reduce memory allocation
            # but also give it the names of the original row order (from the lpmatrix rownames)
            se <- setNames(sapply(1:nrow(X), function(i) {
              X <- X[i, , drop = FALSE]
              sqrt(X %*% Vp %*% t(X)) 
            }), rownames(X)) 
          }
        }
        # ---------------------------------------------------------------
        # Generate predictions without RNS
        # ---------------------------------------------------------------
        
        # Generate the lpmatrix
        X_noRNS <- predict(gmod, newdata = newdata_noRNS, type = "lpmatrix")
        
        # Keep copy for lmpatrix export which has no columns removed
        X_noRNS_full<-X_noRNS 
        
        # Find which bits of the lpmatrix you want to drop (i.e. not predict with) 
        # Based on the columns provided in the function call, any random effects 
        # (determined from the model at the start) AND any RNS bases because these
        # levels are not within the model so predicting without RNS
        if(is.null(to_drop_col)){
          # If nothing is provide in the to_drop_col just list the the random effects
          # and RNS columns
          drop_names_noRNS<-c(grep(RNS_col, colnames(X_noRNS), value=T), random.effects.names)
        } else {
          # If to_drop_col is provided create a vector of the columns in the lpmatrix which
          # refer to that column (collapse "|" allows a vector of multiple columns to drop
          # to be provided in the function call as these will be looked for separately), 
          # any random effects and the RNS 
          drop_names_noRNS<-c(grep(RNS_col, colnames(X_noRNS), value=T), random.effects.names,
                              grep(paste(to_drop_col,collapse="|"), colnames(X_noRNS), value=T))
        }
        
        # Remove any columns you do not want to predict with from the lpmatrix 
        # and model coefficients
        if(!is.null(drop_names_noRNS)){
          # Remove all the lpmatrix you don't want to predict with i.e. the RNS and random effect
          X_noRNS <- X_noRNS[, ! colnames(X_noRNS) %in% drop_names_noRNS]
          # Remove all the coefficients you don't want to predict with i.e. the RNS and random effect
          coefs_noRNS <- coef(gmod)[setdiff(names(coef(gmod)), drop_names_noRNS)]
        }

        # ------------------------------
        # GENERATE THE RESULTS ASKED FOR IN THE FUNCTION CALL 
        # ------------------------------
        if (type == "lpmatrix") {
          cat("RNS lpmatrix generated")
          # return(X) # CANNOT HAVE THIS HERE AS TOO EARLY HAVE TO RETURN 
          # WHEN EVERYTHING HAS BEEN STUCK BACK TOGETHER
        } else if (type == "link") {
          # Just return fit
          #fit <- as.vector(X %*% coefs)
          # Generate the fit but also give it the names of the original 
          # row order (from the lpmatrix rownames)
          fit_noRNS <- setNames(as.vector(X_noRNS %*% coefs_noRNS), rownames(X_noRNS))
          
          if (se.fit) { # if se.fit=FALSE (default) just return fit 
            # If se.fit is true get the standard errors too:
            # Get positions of the columns/rows we don't want
            drop_names_pos_noRNS<-which(names(gmod$coefficients) %in% drop_names_noRNS)
            # Remove columns in Vp (covariance matrix) which are RNS or random effects 
            Vp_noRNS <- gmod$Vp[-drop_names_pos_noRNS, -drop_names_pos_noRNS]
            # Get standard errors row by row to reduce memory allocation
            # but also give it the names of the original row order (from the lpmatrix rownames)
            se_noRNS <- setNames(sapply(1:nrow(X_noRNS), function(i) {
              X_noRNS <- X_noRNS[i, , drop = FALSE]
              sqrt(X_noRNS %*% Vp_noRNS %*% t(X_noRNS)) 
            }), rownames(X_noRNS)) 
          }
        }
        
        # -----------------------------------------------------
        # JOIN OUTPUTS OF THE RNS CATCHMENTS AND NONE RNS CATCHMENTS TOGETHER 
        # -----------------------------------------------------
        # ----
        # ALL FITS
        # ----
        # JOIN THE FITS OF THE RNS CATCHMENTS AND NONE RNS CATCHMENTS TOGETHER 
        all_fit<-c(fit, fit_noRNS)
        
        # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
        # CORRECT SITE
        all_fit<-all_fit[order(factor(names(all_fit), levels = org_order))]
        
        # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
        all_fit<-as.vector(all_fit)
        
        # ----
        # ALL STANDARD ERRORS
        # ----
        if (se.fit) { 
          # JOIN THE SES OF THE RNS CATCHMENTS AND NONE RNS CATCHMENTS TOGETHER 
          all_se<-c(se, se_noRNS)
          
          # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
          # CORRECT SITE
          all_se<-all_se[order(factor(names(all_se), levels = org_order))]
          
          # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
          all_se<-as.vector(all_se)
        }
        # ----
        # ALL LPMATRIX
        # ----
        # JOIN THE SES OF THE RNS CATCHMENTS AND NONE RNS CATCHMENTS TOGETHER 
        all_X<-rbind(X_full, X_noRNS_full)
        
        # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
        # CORRECT SITE
        all_X<-all_X[order(factor(rownames(all_X), levels = org_order)),]
        # CHECK MATRIX
        
        # ------------------------------------------------------
        # NOW NEED TO ADD THE IF STATEMENTS IN FOR WHAT SHOULD BE RETURNED BASED ON FUCNTION CALLS 
        # ------------------------------------------------------
        
        # WRITE OUT THE RESULTS REQUESTED
        if (type == "lpmatrix") {
          return(all_X) 
        } else if (type == "link") {
          if (!se.fit) # if se.fit=FALSE (default) just return fit 
            return(all_fit) 
          return(list(fit = all_fit, se.fit = all_se)) 
        }
        
      }
      
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # GET PREDICTIONS ONLY WITH RIVER NETWORK SMOOTHERS 
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      
      if(nrow(newdata_RNS) > 0 & nrow(newdata_noRNS) == 0) {
        cat("Predictions only with RNS")    

        # Generate the lpmatrix 
        # We know that the RNS bits of this will be wrong and thus we need to update them 
        X <- predict(gmod, newdata = newdata_RNS, type = "lpmatrix")

        # Find appropriate bits of the smooth parts of model (held in gmod$smooth)
        # that are RNS bits using the relevant column that refers to the RNS -in our case 'riversmooth'
        which.smooth <- grep(RNS_col, sapply(gmod$smooth, "[[", "term")) 
        
        # Extract the smoother info for each RNS 
        # This will be a list if the RNS is by a factor (as each are held in a separate list element)
        smooth<-gmod$smooth[which.smooth]
        
        # Get names of all the RNS coefficients
        # The point values (pulled from the rank) relate to the k value provided e.g. 
        # if RNS has K value of 3 (3 df) the RNS will be RNS:FACTOR.1, RNS:FACTOR.2 RNS:FACTOR.3
        coef.names <- unlist(lapply(smooth, function(x) paste0(x$label, ".", 1:x$rank)))
        
        # Find all the RNS covariate names i.e. bits of data used as RNS (
        # in our case Order_Smooth_1 Order_Smooth_2 etc) and reduce the column names 
        # to the number of columns that can be used in the model i.e. to the K value (pulled from the rank)
        # e.g. if RNS has K value of 3 (3 df) the RNS will be Order_Smooth_1, Order_Smooth_2 Order_Smooth_3
        covar.names <-unlist(lapply(smooth, function(x) colnames(x$xt$X)[1:x$rank + ncol(x$xt$X)- x$rank]))
        
        # Replace the elements in the lpmatrix relating to the RNS bases with the appropriate values 
        X[,coef.names] <- as.matrix(newdata_RNS[, covar.names])
        
        # Keep copy for lmpatrix export which has correct RNS columns but no columns removed
        X_full<-X
        
        # Find which bits of the lpmatrix you want to drop (i.e. not predict with) 
        # Based on the columns provided in the function call and any random effects 
        # (determined from the model at the start)
        if(is.null(to_drop_col)){
          # If nothing is provide in the to_drop_col just list the random effects
          drop_names<-random.effects.names
        } else {
          # If to_drop_col is provided create a vector of the columns in the lpmatrix which
          # refer to that column (collapse "|" allows a vector of multiple columns to drop
          # to be provided in the function call as these will be looked for separately) and 
          # any random effects
          drop_names<-c(random.effects.names, grep(paste(to_drop_col,collapse="|"), colnames(X), value=T))
        }
        
        # Remove any columns you do not want to predict with from the lpmatrix 
        # and model coefficients
        if(!is.null(drop_names)){
          # Remove all the lpmatrix you don't want to predict with i.e. the random effects
          X <- X[, ! colnames(X) %in% drop_names]
          # Remove all the coefficients you don't want to predict with i.e. the random effects
          coefs <- coef(gmod)[setdiff(names(coef(gmod)), drop_names)]
        } else{
          # If drop_names is null that means the model has no random effects and no columns which 
          # you do not want to predict with, so you do not need to update the lpmatrix or coefficients 
          coefs <- coef(gmod)
        }
        
        # ------------------------------
        # GENERATE THE RESULTS ASKED FOR IN THE FUNCTION CALL 
        # ------------------------------
        if (type == "lpmatrix") {
          cat("RNS lpmatrix generated") # Don't need to do anything lpmatrix is generated already
        } else if (type == "link") {
          # Generate the fit but also give it the names of the original 
          # row order (from the lpmatrix rownames)
          fit <- setNames(as.vector(X %*% coefs), rownames(X))
          
          if (se.fit) {# if se.fit=FALSE (default) just return fit 
            # If se.fit is true get the standard errors too:
            if(!is.null(drop_names)){
              # Remove any columns you do not want to predict with from the lpmatrix 
              # and model coefficients
              # Get positions of the columns/rows we don't want
              drop_names_pos<-which(names(gmod$coefficients) %in% drop_names)
              # Remove columns in Vp (covariance matrix) that you don't want to predict with
              Vp <- gmod$Vp[-drop_names_pos, -drop_names_pos]
            } else {
              # If drop_names is null that means the model has no random effects and no columns which 
              # you do not want to predict with, so you do not need to update the lpmatrix or coefficients 
              Vp <- gmod$Vp 
            }
            # Get standard errors row by row to reduce memory allocation
            # but also give it the names of the original row order (from the lpmatrix rownames)
            se <- setNames(sapply(1:nrow(X), function(i) {
              X <- X[i, , drop = FALSE]
              sqrt(X %*% Vp %*% t(X)) 
            }), rownames(X)) 
          }
        }
        
        # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
        # CORRECT SITE
        fit<-fit[order(factor(names(fit), levels = org_order))]
        
        # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
        fit<-as.vector(fit)
        
        # ----
        # ALL STANDARD ERRORS
        # ----
        if (se.fit) { 
          # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
          # CORRECT SITE
          se<-se[order(factor(names(se), levels = org_order))]
          
          # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
          se<-as.vector(se)
        }
        
        # ----
        # ALL LPMATRIX
        # ----
        # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
        # CORRECT SITE
        X_full<-X_full[order(factor(rownames(X_full), levels = org_order)),]
        # CHECK MATRIX
        
        # ------------------------------------------------------
        # NOW NEED TO ADD THE IF STATEMENTS IN FOR WHAT SHOULD BE RETURNED BASED ON FUCNTION CALLS 
        # ------------------------------------------------------
        
        # WRITE OUT THE RESULTS REQUESTED
        if (type == "lpmatrix") {
          return(X_full) 
        } else if (type == "link") {
          if (!se.fit) # if se.fit=FALSE (default) just return fit 
            return(fit) 
          return(list(fit = fit, se.fit = se)) 
        }
        
      }
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # GET PREDICTIONS ONLY WITHOUT RIVER NETWORK SMOOTHERS 
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      # ------------------------------------
      
      if(nrow(newdata_RNS) == 0 & nrow(newdata_noRNS) > 0) {
        cat("Predictions only without RNS")  
        
        # Generate the lpmatrix
        X_noRNS <- predict(gmod, newdata = newdata_noRNS, type = "lpmatrix")
        
        # Keep copy for lmpatrix export
        X_noRNS_full<-X_noRNS
        
        # Find which bits of the lpmatrix you want to drop (i.e. not predict with) 
        # Based on the columns provided in the function call, any random effects 
        # (determined from the model at the start) AND any RNS bases because these
        # levels are not within the model so predicting without RNS
        if(is.null(to_drop_col)){
          # If nothing is provide in the to_drop_col just list the the random effects
          # and RNS columns
          drop_names_noRNS<-c(grep(RNS_col, colnames(X_noRNS), value=T), random.effects.names)
        } else {
          # If to_drop_col is provided create a vector of the columns in the lpmatrix which
          # refer to that column (collapse "|" allows a vector of multiple columns to drop
          # to be provided in the function call as these will be looked for separately), 
          # any random effects and the RNS 
          drop_names_noRNS<-c(grep(RNS_col, colnames(X_noRNS), value=T), random.effects.names,
                              grep(paste(to_drop_col,collapse="|"), colnames(X_noRNS), value=T))
        }
        
        # Remove any columns you do not want to predict with from the lpmatrix 
        # and model coefficients
        if(!is.null(drop_names_noRNS)){
          # Remove all the lpmatrix you don't want to predict with i.e. the RNS and random effect
          X_noRNS <- X_noRNS[, ! colnames(X_noRNS) %in% drop_names_noRNS]
          # Remove all the coefficients you don't want to predict with i.e. the RNS and random effect
          coefs_noRNS <- coef(gmod)[setdiff(names(coef(gmod)), drop_names_noRNS)]
        } 
        
        # ------------------------------
        # GENERATE THE RESULTS ASKED FOR IN THE FUNCTION CALL 
        # ------------------------------
        if (type == "lpmatrix") {
          cat("RNS lpmatrix generated")# Don't need to do anything lpmatrix is generated already
        } else if (type == "link") {
          # Generate the fit but also give it the names of the original 
          # row order (from the lpmatrix rownames)
          fit_noRNS <- setNames(as.vector(X_noRNS %*% coefs_noRNS), rownames(X_noRNS))
          
          if (se.fit) { # if se.fit=FALSE (default) just return fit 
            # If se.fit is true get the standard errors too:
            # Remove any columns you do not want to predict with from the lpmatrix 
            # and model coefficients
            # Get positions of the columns/rows we don't want
            drop_names_pos_noRNS<-which(names(gmod$coefficients) %in% drop_names_noRNS)
            # Remove columns in Vp (covariance matrix) which are RNS or random effects 
            Vp_noRNS <- gmod$Vp[-drop_names_pos_noRNS, -drop_names_pos_noRNS]
            # Get standard errors row by row to reduce memory allocation
            # but also give it the names of the original row order (from the lpmatrix rownames)
            se_noRNS <- setNames(sapply(1:nrow(X_noRNS), function(i) {
              X_noRNS <- X_noRNS[i, , drop = FALSE]
              sqrt(X_noRNS %*% Vp_noRNS %*% t(X_noRNS)) 
            }), rownames(X_noRNS)) 
          }
        }
        
        # ----
        # ALL FITS
        # ----
        # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
        # CORRECT SITE
        fit_noRNS<-fit_noRNS[order(factor(names(fit_noRNS), levels = org_order))]
        
        # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
        fit_noRNS<-as.vector(fit_noRNS)
        
        # ----
        # ALL STANDARD ERRORS
        # ----
        if (se.fit) { 
          # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
          # CORRECT SITE
          se_noRNS<-se_noRNS[order(factor(names(se_noRNS), levels = org_order))]
          
          # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
          se_noRNS<-as.vector(se_noRNS)
        }
        
        # ----
        # ALL LPMATRIX
        # ----
        # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
        # CORRECT SITE
        X_noRNS_full<-X_noRNS_full[order(factor(rownames(X_noRNS_full), levels = org_order)),]
        # CHECK MATRIX
        
        # ------------------------------------------------------
        # NOW NEED TO ADD THE IF STATEMENTS IN FOR WHAT SHOULD BE RETURNED BASED ON FUCNTION CALLS 
        # ------------------------------------------------------
        
        # WRITE OUT THE REQUESTED RESULTS
        if (type == "lpmatrix") {
          return(X_noRNS_full) 
        } else if (type == "link") {
          if (!se.fit) # if se.fit=FALSE (default) just return fit 
            return(fit_noRNS) 
          return(list(fit = fit_noRNS, se.fit = se_noRNS)) 
        }
        
      } 
  
    }
    
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # IF MORE THAN ONE RNS FACTOR LEVEL
    # E.G. MULTIPLE CATCHMENTS OR LIFESTAGES
    # SINGLE_RNS = FALSE
    # single_RNS_factor_col = NULL
    # single_RNS_factor = NULL
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    
    # ------------------------------------  
    # Subset to data which you want to predict with and without the RNS 
    # ------------------------------------  
    # Get data with RNS 
    newdata_RNS<-newdata[as.character(newdata[, c(by_factor_col)]) %in% 
                           unlist(lapply(gmod$smooth, "[[", "by.level")),]
    
    # Get data without RNS
    # This could end up being empty if you are not make predictions where you don't have RNS 
    # e.g. only predicting in catchments with data or only predicting for lifestages the model 
    # was fitted to
    newdata_noRNS<-newdata[! as.character(newdata[, c(by_factor_col)]) %in% 
                             unlist(lapply(gmod$smooth, "[[", "by.level")),]
    
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # GET PREDICTIONS WITH & WITHOUT RIVER NETWORK SMOOTHERS 
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    
    if(nrow(newdata_RNS) > 0 & nrow(newdata_noRNS) > 0) {
    cat("Predicting with AND without RNS")
    
    # Generate the lpmatrix 
    # We know that the RNS bits of this will be wrong and thus we need to update them 
    X <- predict(gmod, newdata = newdata_RNS, type = "lpmatrix")

    # Find appropriate bits of the smooth parts of model (held in gmod$smooth)
    # that are RNS bits using the relevant column that refers to the RNS -
    # in our case 'riversmooth'- as RNS has levels everything has to be done with lists now
    which.smooth <- grep(RNS_col, sapply(gmod$smooth, "[[", "term")) 
    
    # Extract the smoother info for each RNS 
    # This will be a list if the RNS is by a factor (as each are held in a separate list element)
    smooth<-gmod$smooth[which.smooth]
    
    # Get names of all the RNS coefficients
    # The point values (pulled from the rank) relate to the k value provided e.g. 
    # if RNS has K value of 3 (3 df) the RNS will be RNS:FACTOR.1, RNS:FACTOR.2 RNS:FACTOR.3
    coef.names <- unlist(lapply(smooth, function(x) paste0(x$label, ".", 1:x$rank)))
    
    # Find all the RNS covariate names i.e. bits of data used as RNS (
    # in our case Order_Smooth_1 Order_Smooth_2 etc) and reduce the column names 
    # to the number of columns that can be used in the model i.e. to the K value (pulled from the rank)
    # e.g. if RNS has K value of 3 (3 df) the RNS will be Order_Smooth_1, Order_Smooth_2 Order_Smooth_3
    covar.names <-lapply(smooth, function(x) colnames(x$xt$X)[1:x$rank + ncol(x$xt$X)- x$rank])
    
    # Give the list elements the name of the Factor of interest 
    names(covar.names)<-lapply(smooth, "[[", "by.level") 
  
    # Get the RNS data for each covariate (RNS base) for each factor level 
    # (fill the list) and save them as a matrix in each list element
    new_data_rns<-lapply(names(covar.names), function(x) {
      id <- newdata_RNS[[by_factor_col]] == x
      as.matrix(newdata_RNS[id, covar.names[[x]]])
    })
    
    # Name the list elements with the coefficient names (rather than the covariate)
    # for easier merging later
    names(new_data_rns)<-lapply(smooth, "[[", "label")
    
    #  Create a matrix to add your new info into
    all_rns_matrix<-matrix(0, nrow=dim(X)[1], ncol=length(coef.names))
    rownames(all_rns_matrix) <- rownames(X)
    colnames(all_rns_matrix) <- coef.names
    
    # Work through each RNS level (each element in the list)
    for(i in 1:length(new_data_rns)){ 
      # Find the position of the rows in matrix with same row names as new_data_rns
      # RNS factor level i
      # Using %in% so don't get warning messages if the vector is not a multiple of all_rns_matrix
      rows<-which(rownames(all_rns_matrix) %in% rownames(new_data_rns[[i]]))
      # create the RNS columns names for subsetting
      cols<-c(paste0(names(new_data_rns[i]),".", 1:ncol(new_data_rns[[i]])))
      # Replace the rows and columns in the matrix of factor level i with the correct rns data
      all_rns_matrix[rows, cols]<-as.matrix(new_data_rns[[i]])
    }  
    
    # Replace the elements in the lpmatrix relating to the RNS bases with the appropriate values  
    X[,coef.names] <- all_rns_matrix
    
    # Keep copy for lmpatrix export which has correct RNS columns but no columns removed
    X_full<-X
    
    # Find which bits of the lpmatrix you want to drop (i.e. not predict with) 
    # Based on the columns provided in the function call and any random effects 
    # (determined from the model at the start) 
    if(is.null(to_drop_col)){
      # If nothing is provided in the to_drop_col just list the random effects
      drop_names<-random.effects.names
    } else {
      # If to_drop_col is provided create a vector of the columns in the lpmatrix which
      # refer to that column (collapse "|" allows a vector of multiple columns to drop
      # to be provided in the function call as these will be looked for separately) and 
      # any random effects
      drop_names<-c(random.effects.names, grep(paste(to_drop_col,collapse="|"), colnames(X), value=T))
    }
    
    # Remove any columns you do not want to predict with from the lpmatrix 
    # and model coefficients
    if(!is.null(drop_names)){
      # Remove all the lpmatrix you don't want to predict with i.e. random effects
      X <- X[, ! colnames(X) %in% drop_names]
      # Remove all the coefficients you don't want to predict with i.e. the RNS and random effect
      coefs <- coef(gmod)[setdiff(names(coef(gmod)), drop_names)]
    } else {
      # If drop_names is null that means the model has no random effects and no columns which 
      # you do not want to predict with, so you do not need to update the lpmatrix or coefficients
      coefs <- coef(gmod)
    }
    
    # ------------------------------
    # GENERATE THE RESULTS ASKED FOR IN THE FUNCTION CALL 
    # ------------------------------
    if (type == "lpmatrix") {
      cat("RNS lpmatrix generated") # Don't need to do anything lpmatrix is generated already
    } else if (type == "link") {
      # Generate the fit but also give it the names of the original 
      # row order (from the lpmatrix rownames)
      fit <- setNames(as.vector(X %*% coefs), rownames(X))
      
      if (se.fit) {# if se.fit=FALSE (default) just return fit 
        # If se.fit is true get the standard errors too:
        if(!is.null(drop_names)){
          # Remove any columns you do not want to predict with from the lpmatrix 
          # and model coefficients
          # Get positions of the columns/rows we don't want
          drop_names_pos<-which(names(gmod$coefficients) %in% drop_names)
          # Remove columns in Vp (covariance matrix) which are RNS or random effects 
          Vp <- gmod$Vp[-drop_names_pos, -drop_names_pos]
        } else {
          # If drop_names is null that means the model has no random effects and no columns which 
          # you do not want to predict with, so you do not need to update the lpmatrix or coefficients 
          Vp <- gmod$Vp 
        }
        # Get standard errors row by row to reduce memory allocation
        # but also give it the names of the original row order (from the lpmatrix rownames)
        se <- setNames(sapply(1:nrow(X), function(i) {
          X <- X[i, , drop = FALSE]
          sqrt(X %*% Vp %*% t(X)) 
        }), rownames(X)) 
      }
    }
    
    # ---------------------------------------------------------------
    # Generate predictions without RNS
    # ---------------------------------------------------------------
    
    # First need to give the RNS smoother factor column a level from the model to generate
    # the lpmatrix (not predicting with the column so can just pick the first in the level)
    # if you don't do this you will get an error saying the value is outsite of the model scope
    newdata_noRNS[by_factor_col]<-unlist(lapply(gmod$smooth, "[[", "by.level"))[1]
    
    # Generate the lpmatrix
    X_noRNS <- predict(gmod, newdata = newdata_noRNS, type = "lpmatrix")
    
    # Keep copy for lmpatrix export
    X_noRNS_full<-X_noRNS
    
    # Find which bits of the lpmatrix you want to drop (i.e. not predict with) 
    # Based on the columns provided in the function call, any random effects 
    # (determined from the model at the start) AND any RNS bases because these
    # levels are not within the model so predicting without RNS
    if(is.null(to_drop_col)){
      # If nothing is provide in the to_drop_col just list the the random effects
      # and RNS columns
      drop_names_noRNS<-c(grep(RNS_col, colnames(X_noRNS), value=T), random.effects.names)
    } else {
      # If to_drop_col is provided create a vector of the columns in the lpmatrix which
      # refer to that column (collapse "|" allows a vector of multiple columns to drop
      # to be provided in the function call as these will be looked for separately), 
      # any random effects and the RNS 
      drop_names_noRNS<-c(grep(RNS_col, colnames(X_noRNS), value=T), random.effects.names,
                          grep(paste(to_drop_col,collapse="|"), colnames(X_noRNS), value=T))
    }
    
    # Remove any columns you do not want to predict with from the lpmatrix 
    # and model coefficients
    if(!is.null(drop_names_noRNS)){
      # Remove all the lpmatrix you don't want to predict with i.e. the RNS and random effect
      X_noRNS <- X_noRNS[, ! colnames(X_noRNS) %in% drop_names_noRNS]
      # Remove all the coefficients you don't want to predict with i.e. the RNS and random effect
      coefs_noRNS <- coef(gmod)[setdiff(names(coef(gmod)), drop_names_noRNS)]
    } 
    
    # ------------------------------
    # GENERATE THE RESULTS ASKED FOR IN THE FUNCTION CALL 
    # ------------------------------
    if (type == "lpmatrix") {
      cat("RNS lpmatrix generated")# Don't need to do anything lpmatrix is generated already
    } else if (type == "link") {
      # Generate the fit but also give it the names of the original 
      # row order (from the lpmatrix rownames)
      fit_noRNS <- setNames(as.vector(X_noRNS %*% coefs_noRNS), rownames(X_noRNS))
      
      if (se.fit) { # if se.fit=FALSE (default) just return fit 
        # If se.fit is true get the standard errors too:
        # Remove any columns you do not want to predict with from the lpmatrix 
        # and model coefficients
        # Get positions of the columns/rows we don't want
        drop_names_pos_noRNS<-which(names(gmod$coefficients) %in% drop_names_noRNS)
        # Remove columns in Vp (covariance matrix) which are RNS or random effects 
        Vp_noRNS <- gmod$Vp[-drop_names_pos_noRNS, -drop_names_pos_noRNS]
        # Get standard errors row by row to reduce memory allocation
        # but also give it the names of the original row order (from the lpmatrix rownames)
        se_noRNS <- setNames(sapply(1:nrow(X_noRNS), function(i) {
          X_noRNS <- X_noRNS[i, , drop = FALSE]
          sqrt(X_noRNS %*% Vp_noRNS %*% t(X_noRNS)) 
        }), rownames(X_noRNS)) 
        
      }
    }
    
    # -----------------------------------------------------
    # JOIN OUTPUTS OF THE RNS CATCHMENTS AND NONE RNS CATCHMENTS TOGETHER 
    # -----------------------------------------------------
    # ----
    # ALL FITS
    # ----
    # JOIN THE FITS OF THE RNS CATCHMENTS AND NONE RNS CATCHMENTS TOGETHER 
    all_fit<-c(fit, fit_noRNS)
    
    # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
    # CORRECT SITE
    all_fit<-all_fit[order(factor(names(all_fit), levels = org_order))]
    
    # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
    all_fit<-as.vector(all_fit)
    
    # ----
    # ALL STANDARD ERRORS
    # ----
    if (se.fit) { 
      # JOIN THE SES OF THE RNS CATCHMENTS AND NONE RNS CATCHMENTS TOGETHER 
      all_se<-c(se, se_noRNS)
      
      # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
      # CORRECT SITE
      all_se<-all_se[order(factor(names(all_se), levels = org_order))]
      
      # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
      all_se<-as.vector(all_se)
    }
    
    # ----
    # ALL LPMATRIX
    # ----
    # JOIN THE SES OF THE RNS CATCHMENTS AND NONE RNS CATCHMENTS TOGETHER 
    all_X<-rbind(X_full, X_noRNS_full)
    
    # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
    # CORRECT SITE
    all_X<-all_X[order(factor(rownames(all_X), levels = org_order)),]
    # CHECK MATRIX
    
    # ------------------------------------------------------
    # NOW NEED TO ADD THE IF STATEMENTS IN FOR WHAT SHOULD BE RETURNED BASED ON FUCNTION CALLS 
    # ------------------------------------------------------
    
    # WRITE OUT THE REQUESTED OUTPUTS
    if (type == "lpmatrix") {
      return(all_X) 
    } else if (type == "link") {
      if (!se.fit) # if se.fit=FALSE (default) just return fit 
        return(all_fit) 
      return(list(fit = all_fit, se.fit = all_se)) 
    }
    
    }
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # GET PREDICTIONS ONLY WITH RIVER NETWORK SMOOTHERS 
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    
    if(nrow(newdata_RNS) > 0 & nrow(newdata_noRNS) == 0) {
      cat("Predicting only with RNS")  
      # ---------------------------------------------------------------
      # Generate predictions with RNS
      # ---------------------------------------------------------------
      
      # Generate the lpmatrix 
      # We know that the RNS bits of this will be wrong and thus we need to update them 
      X <- predict(gmod, newdata = newdata_RNS, type = "lpmatrix")
      
      # Find appropriate bits of the smooth parts of model (held in gmod$smooth)
      # that are RNS bits using the relevant column that refers to the RNS -in our case 'riversmooth'
      which.smooth <- grep(RNS_col, sapply(gmod$smooth, "[[", "term")) 
      
      # Extract the smoother info for each RNS 
      # This will be a list if the RNS is by a factor (as each are held in a separate list element)
      smooth<-gmod$smooth[which.smooth]
      
      # Get names of all the RNS coefficients
      # The point values (pulled from the rank) relate to the k value provided e.g. 
      # if RNS has K value of 3 (3 df) the RNS will be RNS:FACTOR.1, RNS:FACTOR.2 RNS:FACTOR.3
      coef.names <- unlist(lapply(smooth, function(x) paste0(x$label, ".", 1:x$rank)))
      
      # Find all the RNS covariate names i.e. bits of data used as RNS (
      # in our case Order_Smooth_1 Order_Smooth_2 etc) and reduce the column names 
      # to the number of columns that can be used in the model i.e. to the K value (pulled from the rank)
      # e.g. if RNS has K value of 3 (3 df) the RNS will be Order_Smooth_1, Order_Smooth_2 Order_Smooth_3
      covar.names <-lapply(smooth, function(x) colnames(x$xt$X)[1:x$rank + ncol(x$xt$X)- x$rank])
      
      # Give the list elements the name of the Factor of interest 
      names(covar.names)<-lapply(smooth, "[[", "by.level") 
      
      # Get the RNS data for each covariate (RNS base) for each factor level 
      # (fill the list) as save them as a matrix in each list element
      new_data_rns<-lapply(names(covar.names), function(x) {
        id <- newdata_RNS[[by_factor_col]] == x
        as.matrix(newdata_RNS[id, covar.names[[x]]])
      })
      
      # Name the list elements with the coefficient names (rather than the covariate)
      # for easier merging later
      names(new_data_rns)<-lapply(smooth, "[[", "label")
      
      #  Create a matrix to add your new info into
      all_rns_matrix<-matrix(0, nrow=dim(X)[1], ncol=length(coef.names))
      rownames(all_rns_matrix) <- rownames(X)
      colnames(all_rns_matrix) <- coef.names
      
      # Work through each RNS level (each element in the list)
      for(i in 1:length(new_data_rns)){ 
        # Find position of the rows in matrix with same row names as new_data_rns
        # factor level i
        # Using %in% so don't get warning messages if the vector is not a multiple of all_rns_matrix
        rows<-which(rownames(all_rns_matrix) %in% rownames(new_data_rns[[i]]))
        # create the RNS columns names for subsetting
        cols<-c(paste0(names(new_data_rns[i]),".", 1:ncol(new_data_rns[[i]])))
        # Replace the rows and columns in the matrix of factor level i with the correct rns data
        all_rns_matrix[rows, cols]<-as.matrix(new_data_rns[[i]])
      }  
      
      # Replace the elements in the lpmatrix relating to the RNS bases with the appropriate values 
      X[,coef.names] <- all_rns_matrix
      
      # Keep copy for lmpatrix export which has correct RNS columns but no columns removed
      X_full<-X
      
      # Find which bits of the lpmatrix you want to drop (i.e. not predict with) 
      # Based on the columns provided in the function call and any random effects 
      # (determined from the model at the start) 
      if(is.null(to_drop_col)){
        # If nothing is provide in the to_drop_col just list the random effects
        drop_names<-random.effects.names
      } else {
        # If to_drop_col is provided create a vector of the columns in the lpmatrix which
        # refer to that column (collapse "|" allows a vector of multiple columns to drop
        # to be provided in the function call as these will be looked for separately) and 
        # any random effects
        drop_names<-c(random.effects.names, grep(paste(to_drop_col,collapse="|"), colnames(X), value=T))
      }
      
      # Remove any columns you do not want to predict with from the lpmatrix 
      # and model coefficients
      if(!is.null(drop_names)){
        # Remove all the lpmatrix you don't want to predict with i.e. the RNS and random effect
        X <- X[, ! colnames(X) %in% drop_names]
        # Remove all the coefficients you don't want to predict with i.e. the RNS and random effect
        coefs <- coef(gmod)[setdiff(names(coef(gmod)), drop_names)]
      } else {
        # If drop_names is null that means the model has no random effects and no columns which 
        # you do not want to predict with, so you do not need to update the lpmatrix or coefficients 
        coefs <- coef(gmod)
      }
      
      # ------------------------------
      # GENERATE THE RESULTS ASKED FOR IN THE FUNCTION CALL 
      # ------------------------------
      if (type == "lpmatrix") {
        cat("RNS lpmatrix generated") # Don't need to do anything lpmatrix is generated already
      } else if (type == "link") {
        # Generate the fit but also give it the names of the original 
        # row order (from the lpmatrix rownames)
        fit <- setNames(as.vector(X %*% coefs), rownames(X))
        
        if (se.fit) {# if se.fit=FALSE (default) just return fit 
          # If se.fit is true get the standard errors too:
          if(!is.null(drop_names)){
            # Remove any columns you do not want to predict with from the lpmatrix 
            # and model coefficients
            # Get positions of the columns/rows we don't want
            drop_names_pos<-which(names(gmod$coefficients) %in% drop_names)
            # Remove columns in Vp (covariance matrix) which are RNS or random effects 
            Vp <- gmod$Vp[-drop_names_pos, -drop_names_pos]
          } else {
            # If drop_names is null that means the model has no random effects and no columns which 
            # you do not want to predict with, so you do not need to update the lpmatrix or coefficients 
            Vp <- gmod$Vp 
          }
          # Get standard errors row by row to reduce memory allocation
          # but also give it the names of the original row order (from the lpmatrix rownames)
          se <- setNames(sapply(1:nrow(X), function(i) {
            X <- X[i, , drop = FALSE]
            sqrt(X %*% Vp %*% t(X)) 
          }), rownames(X)) 
        }
      }
  
      # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
      # CORRECT SITE
      fit<-fit[order(factor(names(fit), levels = org_order))]
      
      # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
      fit<-as.vector(fit)
      
      # ----
      # ALL STANDARD ERRORS
      # ----
      if (se.fit) { 
        # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
        # CORRECT SITE
        se<-se[order(factor(names(se), levels = org_order))]
        
        # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
        se<-as.vector(se)
      }
      
      # ----
      # ALL LPMATRIX
      # ----
      # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
      # CORRECT SITE
      X_full<-X_full[order(factor(rownames(X_full), levels = org_order)),]
      # CHECK MATRIX
      
      # ------------------------------------------------------
      # NOW NEED TO ADD THE IF STATEMENTS IN FOR WHAT SHOULD BE RETURNED BASED ON FUCNTION CALLS 
      # ------------------------------------------------------
      
      # WRITE OUT REQUESTED OUTPUTS
      if (type == "lpmatrix") {
        return(X_full) 
      } else if (type == "link") {
        if (!se.fit) # if se.fit=FALSE (default) just return fit 
          return(fit) 
        return(list(fit = fit, se.fit = se)) 
      }
      
    }
    
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # GET PREDICTIONS ONLY WITHOUT RIVER NETWORK SMOOTHERS 
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    # ------------------------------------
    
    if(nrow(newdata_RNS) == 0 & nrow(newdata_noRNS) > 0) {
    cat("Predicting only without RNS") 
  
      # First need to give the RNS smoother factor column a level from the model to generate
      # the lpmatrix (not predicting with the column so can just pick the first in the level)
      # if you don't do this you will get an error saying the value is outsite of the model scope
      newdata_noRNS[by_factor_col]<-unlist(lapply(gmod$smooth, "[[", "by.level"))[1]
      
      # Generate the lpmatrix
      X_noRNS <- predict(gmod, newdata = newdata_noRNS, type = "lpmatrix")
      
      # Keep copy for lmpatrix export which has no columns removed
      X_noRNS_full<-X_noRNS
      
      # Find which bits of the lpmatrix you want to drop (i.e. not predict with) 
      # Based on the columns provided in the function call, any random effects 
      # (determined from the model at the start) AND any RNS bases because these
      # levels are not within the model so predicting without RNS
      if(is.null(to_drop_col)){
        # If nothing is provide in the to_drop_col just list the the random effects
        # and RNS columns
        drop_names_noRNS<-c(grep(RNS_col, colnames(X_noRNS), value=T), random.effects.names)
      } else {
        # If to_drop_col is provided create a vector of the columns in the lpmatrix which
        # refer to that column (collapse "|" allows a vector of multiple columns to drop
        # to be provided in the function call as these will be looked for separately), 
        # any random effects and the RNS 
        drop_names_noRNS<-c(grep(RNS_col, colnames(X_noRNS), value=T), random.effects.names,
                            grep(paste(to_drop_col,collapse="|"), colnames(X_noRNS), value=T))
      }
      
      # Remove any columns you do not want to predict with from the lpmatrix 
      # and model coefficients
      if(!is.null(drop_names_noRNS)){
        # Remove all the lpmatrix you don't want to predict with i.e. the RNS and random effect
        X_noRNS <- X_noRNS[, ! colnames(X_noRNS) %in% drop_names_noRNS]
        # Remove all the coefficients you don't want to predict with i.e. the RNS and random effect
        coefs_noRNS <- coef(gmod)[setdiff(names(coef(gmod)), drop_names_noRNS)]
      }
      
      # ------------------------------
      # GENERATE THE RESULTS ASKED FOR IN THE FUNCTION CALL 
      # ------------------------------
      if (type == "lpmatrix") {
        cat("RNS lpmatrix generated") # Don't need to do anything lpmatrix is generated already
      } else if (type == "link") {
        # Generate the fit but also give it the names of the original 
        # row order (from the lpmatrix rownames)
        fit_noRNS <- setNames(as.vector(X_noRNS %*% coefs_noRNS), rownames(X_noRNS))
        
        if (se.fit) { # if se.fit=FALSE (default) just return fit 
          # If se.fit is true get the standard errors too:
          # Remove any columns you do not want to predict with from the lpmatrix 
          # and model coefficients
          # Get positions of the columns/rows we don't want
          drop_names_pos_noRNS<-which(names(gmod$coefficients) %in% drop_names_noRNS)
          # Remove columns in Vp (covariance matrix) which are RNS or random effects 
          Vp_noRNS <- gmod$Vp[-drop_names_pos_noRNS, -drop_names_pos_noRNS]
          # Get standard errors row by row to reduce memory allocation
          # but also give it the names of the original row order (from the lpmatrix rownames)
          se_noRNS <- setNames(sapply(1:nrow(X_noRNS), function(i) {
            X_noRNS <- X_noRNS[i, , drop = FALSE]
            sqrt(X_noRNS %*% Vp_noRNS %*% t(X_noRNS)) 
          }), rownames(X_noRNS)) 
        }
      }
      
      # ----
      # ALL FITS
      # ----
      # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
      # CORRECT SITE
      fit_noRNS<-fit_noRNS[order(factor(names(fit_noRNS), levels = org_order))]
      
      # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
      fit_noRNS<-as.vector(fit_noRNS)
      
      # ----
      # ALL STANDARD ERRORS
      # ----
      if (se.fit) { 
        # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
        # CORRECT SITE
        se_noRNS<-se_noRNS[order(factor(names(se_noRNS), levels = org_order))]
        
        # ONCE IN THE CORRECT ORDER CAN THEN CHANGE BACK TO A VECTOR
        se_noRNS<-as.vector(se_noRNS)
      }
      
      # ----
      # ALL LPMATRIX
      # ----
      # REOORDER THEM TO BE THEIR ORIGINAL COLLUMN ORDER SO THE CORRECT PREDICTIONS GO TO THE 
      # CORRECT SITE
      X_noRNS_full<-X_noRNS_full[order(factor(rownames(X_noRNS_full), levels = org_order)),]
  
      # ------------------------------------------------------
      # NOW NEED TO ADD THE IF STATEMENTS IN FOR WHAT SHOULD BE RETURNED BASED ON FUCNTION CALLS 
      # ------------------------------------------------------
      
      # WRITE OUT REQUESTED OUTPUTS
      if (type == "lpmatrix") {
        return(X_noRNS_full) 
      } else if (type == "link") {
        if (!se.fit) # if se.fit=FALSE (default) just return fit 
          return(fit_noRNS) 
        return(list(fit = fit_noRNS, se.fit = se_noRNS)) 
      }
      
    } 
  }

#save(predict_RNS_gam_bam, file="J:/Schedule_of_Service/FW02G0/Scotland_River_Temperature_Monitoring_Network/National_Tw_Model/R_Code/Multi_catch_RNS_predict_function_dev/predict_RNS_gam_bam.RData")
