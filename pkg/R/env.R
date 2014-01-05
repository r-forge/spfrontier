.SpfrontierEnv = new.env(hash = TRUE)

get.envir = function(){
  res = NULL
  if (exists("estimator.env", envir = .SpfrontierEnv)){
    res = get("estimator.env", envir = .SpfrontierEnv)
  }
  return(res)
}

envir.clone = function(envir){
  return(as.environment(as.list(envir, all.names=TRUE)))
}

envir.init = function(...){
  envir <- get.envir()
  newenv <- new.env()
  if (!is.null(envir)){
    en <- envir.clone(envir)
    assign("par.env", en, envir = newenv)
  }

  argg <- list(...)
  for (a in names(argg)){
    if (!is.null(argg[[a]])) assign(a, argg[[a]], envir=newenv)
  }
  assign("estimator.env", newenv, envir = .SpfrontierEnv)
}

envir.finalise = function(){
  parent.env = envir.get("par.env")
  if(!is.null(parent.env)){
    assign("estimator.env", parent.env, envir = .SpfrontierEnv)
  }else{
    rm("estimator.env",envir = .SpfrontierEnv)
  }
}

envir.get = function(name,...){
  res = NULL
  if (exists(name, envir = get.envir(),...)){
    res = get(name, envir = get.envir(),...)
  }
  return(res)
}

envir.assign = function(name,obj,...){
  assign(name, obj, envir = get.envir(),...)
}

envir.counter = function(name,...){
  suffix = ".count"
  cName = paste(name,suffix,sep="")
  counter = envir.get(cName,...)
  if (!is.null(counter)){
    counter = counter + 1
  }else{
    counter = 1
  }
  assign(cName, counter, envir = get.envir(),...)
  return(counter)
}