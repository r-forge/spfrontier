spfrontier.env <- new.env(hash=TRUE)

spfrontier.env.get = function(name,...){
  res = NULL
  if (exists(name, envir = spfrontier.env,...)){
    res = get(name, envir = spfrontier.env,...)
  }
  return(res)
}

spfrontier.env.clear = function(){
  rm(list = ls(spfrontier.env), pos = spfrontier.env)
}

spfrontier.env.counter = function(name,...){
  suffix = ".count"
  cName = paste(name,suffix,sep="")
  counter = spfrontier.env.get(cName,...)
  if (!is.null(counter)){
    counter = counter + 1
  }else{
    counter = 1
  }
  assign(cName, counter, envir = spfrontier.env,...)
  return(counter)
}