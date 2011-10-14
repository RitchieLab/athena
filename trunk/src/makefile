###### ATHENA Makefile ######

all:
	cd stat_calc; make
	cd data_manage; make
	cd annie; make
	cd galib247; make lib
	cd libGE-0.26; ./configure --with-MITGALIBINCLUDES=../galib247/ga/ --with-MITGALIBLIBS=../galib247/ga/libga.a --enable-GE_ILLIGALSGA=no; make 
	cd athena; make final

clean:
	cd libGE-0.26; make clean
	cd galib247; make clean
	cd annie; make clean
	cd data_manage; make clean
	cd stat_calc; make clean
	cd athena; make clean
