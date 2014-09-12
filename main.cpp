//SDFGen - A simple grid-based signed distance field (level set) generator for triangle meshes.
//Written by Christopher Batty (christopherbatty@yahoo.com, www.cs.columbia.edu/~batty)
//...primarily using code from Robert Bridson's website (www.cs.ubc.ca/~rbridson)
//This code is public domain. Feel free to mess with it, let me know if you like it.

#include "makelevelset3.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <zlib.h>

int main(int argc, char* argv[]) {
  
  if(argc != 4) {
    std::cout << "SDFGen - A utility for converting closed oriented triangle meshes into grid-based signed distance fields.\n";
    std::cout << "\nThe output file format is:";
    std::cout << "<ni> <nj> <nk>\n";
    std::cout << "<origin_x> <origin_y> <origin_z>\n";
    std::cout << "<dx>\n";
    std::cout << "<value_1> <value_2> <value_3> [...]\n\n";
    
    std::cout << "(ni,nj,nk) are the integer dimensions of the resulting distance field.\n";
    std::cout << "(origin_x,origin_y,origin_z) is the 3D position of the grid origin.\n";
    std::cout << "<dx> is the grid spacing.\n\n";
    std::cout << "<value_n> are the signed distance data values, in ascending order of i, then j, then k.\n";

    std::cout << "The output filename will match that of the input, with the OBJ suffix replaced with SDF.\n\n";

    std::cout << "Usage: SDFGen <filename> <dx> <padding>\n\n";
    std::cout << "Where:\n";
    std::cout << "\t<filename> specifies a Wavefront OBJ (text) file representing a *triangle* mesh (no quad or poly meshes allowed). File must use the suffix \".obj\".\n";
    std::cout << "\t<dx> specifies the length of grid cell in the resulting distance field.\n";
    std::cout << "\t<padding> specifies the number of cells worth of padding between the object bound box and the boundary of the distance field grid. Minimum is 1.\n\n";
    
    exit(-1);
  }

  std::string filename(argv[1]);
  if(filename.size() < 5 || filename.substr(filename.size()-4) != std::string(".obj")) {
    std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
    exit(-1);
  }

  std::stringstream arg2(argv[2]);
  float dx;
  arg2 >> dx;
  
  std::stringstream arg3(argv[3]);
  int padding;
  arg3 >> padding;

  if(padding < 1) padding = 1;
  //start with a massive inside out bound box.
  Vec3f min_box(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max()), 
    max_box(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
  
  std::cout << "Reading data.\n";

  std::ifstream infile(argv[1]);
  if(!infile) {
    std::cerr << "Failed to open. Terminating.\n";
    exit(-1);
  }

  int ignored_lines = 0;
  std::string line;
  std::vector<Vec3f> vertList;
  std::vector<Vec3ui> faceList;
  while(!infile.eof()) {
    std::getline(infile, line);
    if(line[0] == 'v' && (line[1] == ' ' || line[1] == '\t')) {
      std::stringstream data(line);
      char c;
      Vec3f point;
      data >> c >> point[0] >> point[1] >> point[2];
      vertList.push_back(point);
      update_minmax(point, min_box, max_box);
    }
    else if(line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
      std::stringstream data(line);
      std::string w;
      int f[4];
      int idx = 0;
      data >> w;
      while (data >> w) {
        std::stringstream wstream(w);
        int n;
        wstream >> n;
        f[idx++] = n;
        if (idx > 4) {
            std::cout << "only tris, quads supported\n";
            exit(1);
        }
      }
      if (idx == 3) {
        faceList.push_back(Vec3ui(f[0]-1,f[1]-1,f[2]-1));
      } else if (idx == 4) {
        faceList.push_back(Vec3ui(f[0]-1,f[1]-1,f[2]-1));
        faceList.push_back(Vec3ui(f[2]-1,f[3]-1,f[0]-1));
      }
    }
    else {
      ++ignored_lines; 
    }
  }
  infile.close();
  
  if(ignored_lines > 0)
    std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";

  std::cout << "Read in " << vertList.size() << " vertices and " << faceList.size() << " faces." << std::endl;

  //Add padding around the box.
  Vec3f unit(1,1,1);
  min_box -= padding*dx*unit;
  max_box += padding*dx*unit;
  Vec3ui sizes = Vec3ui((max_box - min_box)/dx);
  
  std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

  std::cout << "Computing signed distance field.\n";
  Array3f phi_grid;
  make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], phi_grid);

  //Very hackily strip off file suffix.
  std::string outname = filename.substr(0, filename.size()-4) + std::string(".sdf");
  std::cout << "Writing results to: " << outname << "\n";

  int version = 0x1;
  gzFile fp = gzopen(outname.c_str(), "w7");
  gzwrite(fp, &version, sizeof(int));
  gzwrite(fp, &phi_grid.ni, sizeof(int));
  gzwrite(fp, &phi_grid.nj, sizeof(int));
  gzwrite(fp, &phi_grid.nk, sizeof(int));
  gzwrite(fp, &min_box[0], sizeof(float) * 3);
  gzwrite(fp, &dx, sizeof(float));
  gzwrite(fp, &phi_grid.a[0], sizeof(float) * phi_grid.a.size());
  gzclose(fp);
  std::cout << "Processing complete.\n";

return 0;
}
