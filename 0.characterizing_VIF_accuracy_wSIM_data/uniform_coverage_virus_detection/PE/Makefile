DIRS = sim100_PE sim150_PE sim50_PE

all:
	@for i in $(DIRS); do \
	echo "Processing $$i..."; \
	(cd $$i; $(MAKE) all) || exit $$?; done



clean:
	@for i in $(DIRS); do \
	echo "Cleaning $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done
