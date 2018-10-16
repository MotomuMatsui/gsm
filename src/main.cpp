/******************************************\
| Graph Splitting Method v2.1 (2018/10/15) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>
#include <unistd.h>
#include <random>   
#include <functional>

using namespace std;

/// format.cpp (File/Data-handling)
extern void readMAT(ifstream&, double*&, int&);
extern void sc2nwk(int* const&, string&, int const&);
extern void addEP(string const&, string&, unordered_map<string, double>&, int const&, int const&);
extern void addLABEL(string const&, string&, string const&, int const&);

/// gs.cpp (Core functions of GS method)
extern void GS(double* const&, int*&, int const&);
extern void EP(double* const&, unordered_map<string, double>&, function<double()>&, int const&);

/// messages.cpp
extern void print_banner();
extern void print_usage(char*&);

/// transitivity.cpp
extern double transitivity(double* const&, int const&);

int main(int argc, char* argv[]){
  
  /*/ Getopt /*/
  int silence = 0;            // -s
  int ep_num  = 0;            // -e
  int seed    = 0;            // -r

  opterr = 0; // default error messages -> OFF
  int opt;
  regex renum(R"(^[\d\.]+$)"); // -e/-r/-t/-m option requires an integer/flout number
  //while ((opt = getopt(argc, argv, "shlve:r:t:m:")) != -1){
  while ((opt = getopt(argc, argv, "shve:r:")) != -1){
    if(opt == 'e'){ // OK! (./gs -e 100 IN.fst)
      if(regex_match(optarg, renum)){
	ep_num = atoi(optarg);
      }
      else{ // NG! (./gs -e hundred IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'r'){ // OK! (./gs -r 12345 IN.fst)
      if(regex_match(optarg, renum)){
	seed = atoi(optarg);
      }
      else{ // NG! (./gs -r one_two_three IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -r requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'h'){ // HELP message (./gs -h)
      /*PRINT*/ print_banner();
      /*PRINT*/ print_usage(argv[0]);
      return 0;
    }
    else if(opt == 'v'){ // Version (./gs -v)
      /*PRINT*/ print_banner();
      return 0;
    }
    else if(opt == 's'){ // SILENT mode (./gs -s -e 100 IN.fst)
      silence = 1;
    }
    else if (opt == '?'){
      if(optopt == 'e'){ // NG! (./gs IN.fst -e)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
      else if(optopt == 'r'){ // NG! (./gs IN.fst -r)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -r requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
      else{ // NG! (./gs -Z)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << argv[0] << ": invalid option\n" <<  endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
  }

  /*/ Input file name /*/
  string input = "";
  if(optind < argc){ // OK! (./gs -e 100 IN.fst)
    input = argv[optind];
  }
  else{ // NG! (./gs -e 100)
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << argv[0] << " requires an input file (matrix).\n" << endl;
    /*PRINT*/ print_usage(argv[0]);
    return -1;
  }
  
  /*/ Variables /*/
  double* W;     // Sequence similarity matrix
  int size;      // Row size of W (W is a synmetry matrix)
  int* gs;       // Result of GS method
  string newick; // GS tree (without EP values)

  /*/File I/O/*/
  auto original_file = string(input);
  ifstream ifs1(original_file); // Matrix file (original)
  if(ifs1.fail()){
    /*PRINT*/ cerr << "\nCannot access " << original_file << "!" << endl;
    return -1;
  }

  /*/Parsing matrix file/*/
  readMAT(ifs1, W, size); 
  // ifs1: INPUT (original matrix file)
  // W:    OUTPUT (matrix)
  // size: # of sequence = row size of sequence similarity matrix

  auto transitivity_score = transitivity(W, size);
  
  /*/ Parameters /*/  
  if(!silence){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "Settings:" << endl;
    /*PRINT*/ cerr << "-Input" << endl;
    /*PRINT*/ cerr << "  File = " << input << endl;
    /*PRINT*/ cerr << "  # of sequences = " << size << endl << endl;
    /*PRINT*/ cerr << "  Transitivity = "<< transitivity_score << endl << endl;

    /*PRINT*/ cerr << "-EP method" << endl;
    if(seed>0){
      /*PRINT*/ cerr << "  Random seed = " << seed << endl;
    }
    else{
      /*PRINT*/ cerr << "  Random seed = " << "a random number (default)" << endl;
    }
    /*PRINT*/ cerr << "  # of iterations = " << ep_num << endl << endl;

    /*PRINT*/ cerr << "Progress:" << endl;
  }

  /*/ GS method (stepwise spectral clustering) /*/
  /*PRINT*/ if(!silence) cerr << "-GS method\n" << "  executing...\r" << flush;
  GS(W, gs, size);
    // W: INPUT (sequence similarity matrix)
    // gs: OUTPUT (result of stepwise spectral clustering)

  /*PRINT*/ if(!silence) cerr << "  done.         " << endl << endl;

  /*/ Generating GS tree Newick based on the spectral clustering /*/
  sc2nwk(gs, newick, size);
    // gs: INPUT (result of stepwise spectral clustering)
    // newick: OUTPUT (GS tree [newick format])

  /*/ EP method /*/
  if(ep_num>0){
    unordered_map<string, double> ep;
    string newick_EP; // GS+EP tree

    // Random number generator (Uniform distribution->Mersenne Twister)
    function<double()> R;
    uniform_real_distribution<double> urd(0,1);    // uniform distributed random number

    if(seed>0){
      mt19937 mt(static_cast<unsigned int>(seed)); // mersenne twister
      R = bind(urd, ref(mt));                      // random number generator    
    }
    else{
      random_device rd;                            // random seed
      mt19937 mt(rd());                            // mersenne twister
      R = bind(urd, ref(mt));                      // random number generator        
    }    

    /*PRINT*/ if(!silence) cerr << "-EP method" << endl;

    for(int n=1; n<=ep_num; n++){
      /*PRINT*/ if(!silence) cerr << "  " << n << "/" << ep_num << " iterations" << "\r"<< flush;
      EP(W, ep, R, size);
        // W: INPUT (sequence similarity matrix)
        // ep: OUTPUT (result of Edge Perturbation method)
        // R: random number generator
    }
    
    /*PRINT*/ if(!silence) cerr << "\n  done." << endl << endl;
    /*PRINT*/ if(!silence) cerr << "------------------------------------------\n" << endl;

    addEP(newick, newick_EP, ep, ep_num, size);
      // newick: INPUT (GS tree [newick format])
      // newick_EP: OUTPUT (GS+EP tree [newick format])
      // ep: INPUT (result of Edge Perturbation method)
      // ep_num: INPUT (# of Edge Perturbation method)

    /*/ GS tree WITH EP values ->STDOUT /*/
    cout << newick_EP << endl;
  }
  else{ // skip the EP method
    /*PRINT*/ if(!silence) cerr << "------------------------------------------\n" << endl;

    /*/ GS tree WITHOUT EP values ->STDOUT /*/
    cout << newick << endl;
  }

  delete[] W;
  delete[] gs;

  return 0;
}
