logger.levels = c("quiet","warn","info","debug")

logger.log = function(message, level = "debug", obj = NULL, caller = sys.call(-1)){
  lvl = spfrontier.env.get("logging.level")
  is.log = !is.null(lvl) && (which(logger.levels == lvl) >=  which(logger.levels == level))
  dt = paste("[",format(Sys.time(), "%c"),"]",sep="")
  if (is.log){
    print(paste(dt, toupper(level), caller[[1]], message))
    if (!is.null(obj)){
      print(obj)
    }
  }
}

logger.warn = function(message, obj = NULL,caller = sys.call(-1)){
  logger.log(message,"warn", obj = obj,caller)
}

logger.info = function(message, obj = NULL, caller = sys.call(-1)){
  logger.log(message,"info", obj = obj,caller)
}

logger.debug = function(message, obj = NULL,caller = sys.call(-1)){
  logger.log(message,"debug", obj = obj,caller)
}