
cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools
module add samtools/1.14

# repair is used to put mates next to each other; 
# it's part of the subread package (http://subread.sourceforge.net/); 
# if there's no module for subread, then you can do a local installation through sourceforge or conda


ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/bin/subread-2.0.3-Linux-x86_64 subread
subread/bin/utilities/repair 

