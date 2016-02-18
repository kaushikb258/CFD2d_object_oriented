all:
	g++ -g main.cpp -o cfd2d.exe
clean:
	rm *.exe *~ output_file *.vtk
