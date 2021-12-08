Getriskset <- function(cdata, surv.formula) {
  
  survname <- all.vars(surv.formula)
  status <- as.vector(cdata[, survname[2]])
  
  subcdata <- cdata[, survname[1:2]]
  subcdata <- as.matrix(subcdata)
  
  if (prod(c(0, 1, 2) %in% unique(status))) {
    
    riskset <- GetrisksetC(subcdata)
    tablerisk1 <- riskset$H01
    tablerisk2 <- riskset$H02
    
    colnames(tablerisk1) <- c("survtime", "d", "hazard")
    colnames(tablerisk2) <- c("survtime", "d", "hazard")
    a <- list(tablerisk1, tablerisk2)
    names(a) <- c("tablerisk1", "tablerisk2")
    return(a)
  } else {
    
    riskset <- GetrisksetCSF(subcdata)
    tablerisk1 <- riskset$H01
    
    colnames(tablerisk1) <- c("survtime", "d", "hazard")
    a <- list(tablerisk1)
    names(a) <- c("tablerisk1")
    return(a)
    
  }
  
}