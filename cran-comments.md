There is one note upon the R CMD check. It is as follows:

  "* checking for unstated dependencies in vignettes ... NOTE
'::' or ':::' imports not declared from:
  ‘BiocManager’ ‘EBImage’ ‘readr’
'library' or 'require' calls not declared from:
  ‘BiocManager’ ‘tidyverse’"

This is due to the use of packages, ‘BiocManager’ ‘EBImage’, ‘readr’, and ‘tidyverse’, in the vignette. These functions are not independently called as part of the package, i.e., no function in the SIPmg R package utilizes directly these packages listed above. 

Since these packages listed above are only part of the vignette, the CRAN generated NOTE can be ignored.
