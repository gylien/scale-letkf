.PHONY:	all clean

all:
	$(MAKE) -C common
	$(MAKE) -C .common_for_scale
	$(MAKE) -C ensmodel
	$(MAKE) -C obs
	$(MAKE) -C letkf

clean:
	$(MAKE) -C common            clean
	$(MAKE) -C .common_for_scale clean
	$(MAKE) -C ensmodel          clean
	$(MAKE) -C obs               clean
	$(MAKE) -C letkf             clean
