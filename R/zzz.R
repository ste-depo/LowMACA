	
.onLoad <- function(libname, pkgname){
    #packageStartupMessage('Welcome to the wonderful World of LowMACA', domain = NULL, appendLF = TRUE)
    	#Check for clustalomega installation and version	
    .ClustalChecks(ClustalCommand="clustalo")
    	#Check for perl modules dependencies
	.PerlModuleChecks(stop=FALSE , perl="perl")

}
