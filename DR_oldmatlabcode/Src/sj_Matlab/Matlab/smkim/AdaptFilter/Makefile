
DIRS= aux stats model gen filter main newvis

all:
	make -C aux
	make -C model
	make -C filter
	make -C gen
	make -C stats
	make -C vis

depend:
	make -C aux depend
	make -C stats depend
	make -C model depend
	make -C filter depend
	make -C gen depend
	make -C vis depend

clean:
	make -C aux clean
	make -C stats clean
	make -C model clean
	make -C filter clean
	make -C gen clean
	make -C vis clean
	make -C tests clean

include makeconfig
