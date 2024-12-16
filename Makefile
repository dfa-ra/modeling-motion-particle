run:
	cmake -S . -B build
	make -C ./build
	eval ./build/projectC
	python3 ./drawer/draw.py