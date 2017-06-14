all:
	$(MAKE) -C ./src/convert all
	cp ./src/convert/convert .
	$(MAKE) -C ./src/rect_trian_mesh
	cp ./src/rect_trian_mesh/rect_trian_mesh .
clean:
	$(MAKE) -C ./src/convert clean
	$(MAKE) -C ./src/rect_trian_mesh clean
