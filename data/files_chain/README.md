# files_chain

```
#!/bin/bash

#  Call script from the repo's home directory, 2021_kga0_4dn-mouse-cross

#  Script will download chain files for the following genome coordinate
#+ conversions:
#+ - CAST-EiJ to mm10
#+ - mm10 to CAST-EiJ
#+ - 129S1-SvImJ to mm10
#+ - mm10 to 129S1-SvImJ
#+ - CAROLI-EiJ to mm10
#+ - mm10 to CAROLI-EiJ
#+ - SPRET-EiJ to mm10
#+ - mm10 to SPRET-EiJ
#+ 
#+ Default settings are for files to be saved in 'data/files_chain' and for
#+ chains to be renamed, e.g., 'GCA_001624445.1ToMm10.over.chain.gz' is renamed
#+ 'CAST-EiJ-to-mm10.over.chain.gz'

#  Call with default settings
bash bin/get-liftOver-chains.sh

#  Run time: 5 sec

#  #DONTRUN Call with user-specified settings (FALSE: don't rename chain files)
#+          $ bash bin/get-liftOver-chains.sh path/to/store/files FALSE
```
