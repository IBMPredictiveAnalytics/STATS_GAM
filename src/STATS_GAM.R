#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2016
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "IBM SPSS, JKP"
# version__ = "1.0.0"

# History
# 31-jul-2016 Original Version


gtxt <- function(...) {
    return(gettext(...,domain="STATS_GAM"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_GAM"))
}

missingmap= list("listwise"=na.exclude, "fail"=na.fail)
### MAIN ROUTINE ###

dogam = function(dependent, variables0=NULL,id=NULL, 
    errordist="gaussian", link="identity",
    variables1=NULL, smoother1=NULL, parm11=NULL, parm12=NULL,
    variables2=NULL, smoother2=NULL, parm21=NULL, parm22=NULL,
    variables3=NULL, smoother3=NULL, parm31=NULL, parm32=NULL,
    variables4=NULL, smoother4=NULL, parm41=NULL, parm42=NULL,
    variables5=NULL, smoother5=NULL, parm51=NULL, parm52=NULL,
    missingaction="listwise", offset=NULL, plotit=TRUE,
    epsilon=1e-07, bfepsilon=1e-07, maxit=30, bfmaxit=30, iterations=FALSE,
    dataset=NULL, ptype="link", stderr=TRUE
    ) {
  # Estimate GAM

  setuplocalization("STATS_GAM")
  
  # A warnings proc name is associated with the regular output
  # (and the same omsid), because warnings/errors may appear in
  # a separate procedure block following the regular output
  procname=gtxt("Generalized Additive Model")
  warningsprocname = gtxt("Generalized Additive Model: Warnings")
  omsid="STATSGAM"
  warns = Warn(procname=warningsprocname,omsid=omsid)

  tryCatch(library(gam), error=function(e){
      warns$warn(gtxtf("The R %s package is required but could not be loaded.", "gam"),dostop=TRUE)
      }
  )

  if (!is.null(spssdata.GetSplitVariableNames())) {
      warns$warn(
          gtxt("Split variables are not honored by this procedure"),
          dostop=FALSE)
  }

  checkdataset(dataset, warns)
  errandlink = makeerrlink(errordist, link, warns)
  alldata = c(dependent, variables0, variables1, variables2, variables3, variables4, variables5, offset)
  allargs = as.list(environment())
  dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE, factorMode="labels", row.label=id)
  dtalen = nrow(dta)
  dta = dta[complete.cases(dta),]
  allargs$missingcount = dtalen - nrow(dta)
  frml = as.formula(buildfrml(allargs, warns))
  res = tryCatch(gam(frml, family=errandlink, data=dta, na.action=missingmap[[missingaction]],
    offset=offset,
    epsilong=epsilon, bf.epsilon=bfepsilon, maxit=maxit, bf.maxit=bfmaxit, trace=iterations),
    error=function(e) {warns$warn(e$message, dostop=TRUE)}
  )
  displayresults(allargs, res, warns)
  makedataset(dataset, res, dta, id, ptype, stderr, warns)

  warns$display()
}

checkdataset = function(dataset, warns) {
  # Ensure that the active dataset is named
  # dataset is the name for a new dataset - must not already exist

  if (is.null(dataset)) {
    return()
  }
    alldatasets = spssdata.GetDataSetList()
  if ("*" %in% alldatasets) {
    warns$warn(gtxt("The active dataset must have a name in order to use this procedure"),
               dostop=TRUE)
  }
  if (dataset %in% alldatasets) {
    warns$warn(gtxt("The output dataset name must not already be in use"),
               dostop=TRUE)
  }
  
}
linkmap = list(identity="identity", log="log", inverse="inverse", 
  logit="logit", probit="probit", cauchit="cauchit", 
  cloglog="cloglog", sqrt="sqrt", mu2inv="1/mu^2")
errlinkmap=list(gaussian=c("identity", "log", "inverse"), binomial=c("logit", "probit", "cauchit", "log", "cloglog"),
  Gamma=c("inverse", "identity", "log"), poisson=c("log", "identity", "sqrt"),
  inverse.gaussian=c("m2inv", "inverse", "identity", "log"),
  quasi=c("logit", "probit", "cloglog", "identity", "inverse", "log", "m2inv", "sqrt")
)
makeerrlink = function(errordist, link, warns) {
  # return a family object with the error distribution
  # and link.
  # validation is via the family function but
  # might not be caught until used
  # errordist is the error distribution name, which cannot be NULL
  # link is the link name and can be NULL, in which case
  # the distribution default is used
  
  
  if (errordist == "gamma") {
    errordist = "Gamma"
  }
  # check that a specified link is valid for the specified error distribution
  if (!is.null(link) && substr(errordist, 1, 5) !="quasi") {
    if (!link %in% errlinkmap[[errordist]] ) {
      warns$warn(gtxtf("The specified link is not valid for the error distribution: %s", link), dostop=TRUE)
    }
  }
  link = linkmap[[link]]
  ff=tryCatch(do.call(errordist, list(link=link)),
    error=function(e) {warns$warn(e$message, dostop=TRUE)}
  )
  return(ff)
}

makedataset = function(dsname, res, dta, id, ptype, stderr, warns) {
  # create dataset for predicted values
  # dsname is the name of the dataset to create
  # res is the gam output
  # dta is the input data
  # id is the name of the id variable, which might be NULL
  # stderr indicates whether to include prediction standard errors
  
  if (is.null(dsname)) {
    return()
  }

  predvalues = predict(res, type=ptype, se.fit=stderr)
  # compact the names
  if (ptype == "terms") {
    dn = dimnames(predvalues$fit)[[2]]
    dn = gsub(" ", "", dn)
    dn = gsub("=", " EQ ", dn)
    dimnames(predvalues$fit)[[2]] = dn
  }
  if (!stderr) {
    preddf = data.frame(row.names=row.names(dta), predvalues$fit)
  } else {
      if (ptype == "terms") {
        dn = dimnames(predvalues$se.fit)[[2]]
        dn = gsub(" ", "", dn)
        dn = gsub("=", " EQ ", dn)
        dimnames(predvalues$se.fit)[[2]] = dn
      }
    preddf = data.frame(row.names=row.names(dta), predvalues$fit, predvalues$se.fit)
}

  # make dictionary for new dataset
  dict = list()
  nn = names(preddf)
  idlen = max(nchar(row.names(preddf), type="bytes"))
  if (!is.null(id)) {
    idname = id
  } else {
    id = "ID"
    for (i in 1:9999) {
      idname = paste(id, i, sep="")
      if (!idname %in% toupper(nn)) {
        break
      }
    }
  }

  dict[[1]] = c(idname, "", idlen, paste("A", idlen, sep=""), "nominal")
  
  for (v in 1:length(nn)) {
    dict[[v+1]] = c(nn[[v]], "", 0, "F8.2", "scale")
  }
  
  dict = spssdictionary.CreateSPSSDictionary(dict)
  spssdictionary.SetDictionaryToSPSS(dsname, dict)
  tryCatch(spssdata.SetDataToSPSS(dsname, data.frame(row.names(preddf), preddf)),
           error=function(e) {warns$warn(e$message, dostop=TRUE)}
  )
  spssdictionary.EndDataStep()
  
}


smootherpatterns = list(s="s(%s, df=%s)", lo="lo(%s, span=%s, degree=%s)", 
  bs="bs(%s, df=%s, degree=%s)", poly="poly(%s, degree=%s)", ns="ns(%s, df=%s)")
numparam = list(s=1, lo=2, bs=2, poly=1, ns=1)

buildfrml = function(alla, warns) {
  # build formula from variables, smoother, parm1, and parm2 sets of arguments
  # return as a string
  # variables is a simple list.  The others require a smoother and one or two parms
  
  vars = list(alla$variables1, alla$variables2, alla$variables3, alla$variables4, alla$variables5)
  smoothers = list(alla$smoother1, alla$smoother2, alla$smoother3, alla$smoother4, alla$smoother5)
  parms1 = list(alla$parm11, alla$parm21, alla$parm31, alla$parm41, alla$parm51)
  parms2 = list(alla$parm12, alla$parm22, alla$parm32, alla$parm42, alla$parm52)

  # start of formula and any linear variables
  f = paste(alla$dependent, "~")

  if (!is.null(alla$variables0)) {
    variables = paste(alla$variables0, collapse="+")
    f = paste(f, variables, collapse="")
  }
  # append independent variable specifications, if any
  i = 0

  for (varset in vars) {
    i = i + 1
    if (is.null(varset)) {
      next
    }
    if (is.null(smoothers[[i]]) || is.null(parms1[[i]])) {
      warns$warn(gtxtf("A smoother or smoother parameter was not specified for a variable set: %s.", i), dostop=TRUE)
    }
    pat = tryCatch(smootherpatterns[[smoothers[[i]]]], 
        error=function(e) {warns$warn(gtxtf("Invalid smoother: %s", smoothers[[i]]), dostop=TRUE)})
    # apply smoother to each variable in the current list
    for (v in varset) {
      if (numparam[[smoothers[[i]]]] == 1) {
        patv = sprintf(pat, v, parms1[[i]])
      } else {
        if (is.null(parms2[[i]])) {
          warns$warn(gtxtf("The second parameter for a smoother was not specified: %s", smoothers[[i]]), dostop=TRUE)
        }
        patv = sprintf(pat, v, parms1[[i]], parms2[[i]])
      }
      if (!substr(f, nchar(f), nchar(f)) %in% list("~", "+")) {
        f = paste(f, "+", collapse="")
      }
      f = paste(f, patv, collapse="+")
    }
  }
  if (substr(f, nchar(f), nchar(f)) == "~") {
    warns$warn(gtxt("No independent variables were specified"), dostop=TRUE)
  }
  return(f)
}

displayresults = function(allargs, results, warns) {
    # display results
    # allargs is the parameter set
    

    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    
    # summary results
    # input specifications
    # although groups can be specified (cengroup), separate results are not
    # produced.
    lbls = c(gtxt("Dependent Variable"),
             gtxt("Error Distribution"),
             gtxt("Link"),
             gtxt("Dispersion"),
             gtxt("Offset"),
             gtxt("Cases deleted for missing data"),
             gtxt("Residual D. F."),
             gtxt("Null Model Residual D. F."),
             gtxt("Deviance"),
             gtxt("Null Deviance"),
             gtxt("AIC"),
             gtxt("BIC"),
             gtxt("Iter"),
             gtxt("Convergence Threshold for Local Scoring Iterations"),
             gtxt("Convergence Threshold for Backfitting Iterations"),
             gtxt("Maximum Number of Local Scoring Iterations"),
             gtxt("Maximum Number of Backfitting Iterations"),
             gtxt("Output Dataset")
    )

    vals = c(
            allargs$dependent,
            results$family$family,
            results$family$link,
            round(summary(results)$dispersion[[1]],4),
            ifelse(is.null(allargs$offset), gtxt("--NA--"), allargs$offset),
            allargs$missingcount,
            round(results$df.residual,4),
            round(results$df.null, 4),
            round(results$deviance, 4),
            round(results$null.deviance, 4),
            round(AIC(results), 4),
            round(BIC(results), 4),
            results$iter,
            allargs$epsilon,
            allargs$bfepsilon,
            allargs$maxit,
            allargs$bfmaxit,
            ifelse(is.null(allargs$dataset), gtxt("--NA--"), allargs$dataset)
      
    )
    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Summary"),
        collabels=c(gtxt("Summary")), templateName="GAMSUMMARY", outline=gtxt("Summary"),
        caption = gtxtf("Computations done by R package gam, version: %s", packageVersion("gam"))
    )

    a = summary(results)$parametric.anova
    if (!is.null(a)) {
      spsspivottable.Display(data.frame(a), title=attr(a, "heading"),
        templateName="PARAMETRICANOVA",
        collabels=list(gtxt("D. F"), gtxt("Sum of Squares"), gtxt("Mean Square"),
        gtxt("F"), gtxt("Significance"))
    )
    }
    # nonparametric anova
    a = anova(results)
    if (!is.null(a)) {
      spsspivottable.Display(a, title=attr(a, "heading"),
        templateName="NONPARAMETRICANOVA", 
        collabels=list(gtxt("Nonpar D. F."), gtxt("Nonpar F"), gtxt("Significance"))
      )
    }

    spsspivottable.Display(data.frame(coefficients(results)), title=gtxt("Coefficients"),
      collabels=gtxt("Coefficients"),
      templateName="GAMCOEFFICIENTS", outline=gtxt("Coefficients"))
    
    if (allargs$plotit) {
      plot(results, se=allargs$stderr, lwd=2, main=gtxt("Terms"))
    }

    spsspkg.EndProcedure()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}



Run = function(args) {
    #Execute the STATS GAM command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("DEPENDENT", subc="", ktype="existingvarlist",
            var="dependent", islist=FALSE),
        spsspkg.Template("VARIABLES", subc="", ktype="existingvarlist", var="variables0", islist=TRUE),
        spsspkg.Template("ID", subc="", ktype="existingvarlist", var="id"),
        spsspkg.Template("ERRORDIST", subc="", ktype="str", var="errordist",
          vallist=list("binomial", "gaussian", "gamma", "inversegaussian",
            "poisson", "quasi", "quasibinomial", "quasipoisson")),
        spsspkg.Template("LINK", subc="", ktype="str", var="link",
          vallist=list("identity", "log", "inverse", "logit", "probit",
          "cauchit", "cloglog", "sqrt", "mu2inv")),
        spsspkg.Template("OFFSET", subc="", ktype="existingvarlist", var="offset"),
        
        spsspkg.Template("VARIABLES", subc="NONLINEAR1", ktype="existingvarlist", var="variables1", islist=TRUE),
        spsspkg.Template("SMOOTHER", subc="NONLINEAR1", ktype="str", var="smoother1",
            vallist=list("s", "lo", "bs", "ns", "poly")),
        spsspkg.Template("PARM1", subc="NONLINEAR1", ktype="float", var="parm11"),
        spsspkg.Template("PARM2", subc="NONLINEAR1", ktype="float", var="parm12"), 
        
        spsspkg.Template("VARIABLES", subc="NONLINEAR2", ktype="existingvarlist", var="variables2", islist=TRUE),
        spsspkg.Template("SMOOTHER", subc="NONLINEAR2", ktype="str", var="smoother2",
                         vallist=list("s", "lo", "bs", "ns", "poly")),
        spsspkg.Template("PARM1", subc="NONLINEAR2", ktype="float", var="parm21"),
        spsspkg.Template("PARM2", subc="NONLINEAR2", ktype="float", var="parm22"), 
        
        spsspkg.Template("VARIABLES", subc="NONLINEAR3", ktype="existingvarlist", var="variables3", islist=TRUE),
        spsspkg.Template("SMOOTHER", subc="NONLINEAR3", ktype="str", var="smoother3",
                         vallist=list("s", "lo", "bs", "ns", "poly")),
        spsspkg.Template("PARM1", subc="NONLINEAR3", ktype="float", var="parm31"),
        spsspkg.Template("PARM2", subc="NONLINEAR3", ktype="float", var="parm32"), 
        
        spsspkg.Template("VARIABLES", subc="NONLINEAR4", ktype="existingvarlist", var="variables4", islist=TRUE),
        spsspkg.Template("SMOOTHER", subc="NONLINEAR4", ktype="str", var="smoother4",
                         vallist=list("s", "lo", "bs", "ns", "poly")),
        spsspkg.Template("PARM1", subc="NONLINEAR4", ktype="float", var="parm41"),
        spsspkg.Template("PARM2", subc="NONLINEAR4", ktype="float", var="parm42"), 
        
        spsspkg.Template("VARIABLES", subc="NONLINEAR5", ktype="existingvarlist", var="variables5", islist=TRUE),
        spsspkg.Template("SMOOTHER", subc="NONLINEAR5", ktype="str", var="smoother5",
                         vallist=list("s", "lo", "bs", "ns", "poly")),
        spsspkg.Template("PARM1", subc="NONLINEAR5", ktype="float", var="parm51"),
        spsspkg.Template("PARM2", subc="NONLINEAR5", ktype="float", var="parm52"),
        
        spsspkg.Template("EPSILON", subc="OPTIONS", ktype="float", var="epsilon"),
        spsspkg.Template("BFEPSILON", subc="OPTIONS", ktype="float", var="bfepsilon"),
        spsspkg.Template("MAXIT", subc="OPTIONS", ktype="float", var="maxit"),
        spsspkg.Template("BFMAXIT", subc="OPTIONS", ktype="float", var="bfmaxit"),
        spsspkg.Template("SHOWITERATIONS", subc="OPTIONS", ktype="bool", var="iterations"),
        spsspkg.Template("PLOT", subc="OPTIONS", ktype="bool", var="plotit"),
        
        spsspkg.Template("ACTION", subc="MISSING", ktype="str", var="missingaction",
          vallist=list("listwise", "fail")),
        
        spsspkg.Template("DATASET", subc="SAVE", ktype="varname", var="dataset"),
        spsspkg.Template("TYPE", subc="SAVE", ktype="str", var="ptype",
          vallist=list("link", "response", "terms")),
        spsspkg.Template("STDERR", subc="SAVE", ktype="bool", var="stderr")
        
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "dogam")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
