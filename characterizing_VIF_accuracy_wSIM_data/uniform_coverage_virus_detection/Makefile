DIRS = PE SE

test:
        @for i in $(DIRS); do \
        echo "Processing $$i..."; \
        (cd $$i; $(MAKE) test) || exit $$?; done



clean:
        @for i in $(DIRS); do \
        echo "Cleaning $$i..."; \
        (cd $$i; $(MAKE) clean) || exit $$?; done
