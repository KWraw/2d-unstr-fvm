
COMPONENTS = 2d-unstr-fvm-adv mesh-utilities manual
.Phony: $(COMPONENTS)$ clean


all: $(COMPONENTS)

2d-unstr-fvm-adv:
	$(MAKE) -C ./src

mesh-utilities:
	$(MAKE) -C ./mesh

manual:
	$(MAKE) -C ./doc

clean:
	$(MAKE) clean -C ./src
	$(MAKE) clean -C ./mesh
	$(MAKE) clean -C ./doc
