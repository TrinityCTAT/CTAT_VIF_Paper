DIRS = sim100_SE sim150_SE sim50_SE

all:
	@for i in $(DIRS); do \
	echo "Processing $$i..."; \
	(cd $$i; $(MAKE) all) || exit $$?; done



clean:
	@for i in $(DIRS); do \
	echo "Cleaning $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done
