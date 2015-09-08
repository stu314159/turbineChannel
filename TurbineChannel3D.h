#ifndef TURBINECHANNEL_H
#define TURBINECHANNEL_H

#include <mpi.h>
#include <string>
#include <set>
#include <map>

#define HALO 1

using namespace std;

class TurbineChannel3D{
public:

	// MPI data
	const int rank,size;
	MPI_Status stat;
	MPI_Request rq_in1,rq_in2;
	MPI_Request rq_out1,rq_out2;
	int tag_d, tag_u;

	//input files
        static const string params_file = "params.lbm";
        static const string inl_file = "inl.lbm";
        static const string onl_file = "onl.lbm";
        static const string snl_file = "snl.lbm";

	//input data
	int LatticeType;
	int Num_ts;
	int ts_rep_freq;
        int Warmup_ts;
	int plot_freq;
        float Cs;
	
	float rho_lbm;
	float umax_lbm;
	float omega;
	int Nx;
	int Ny;
	int Nz;

	// lattice data
	static const int numSpd=15;
	const float * ex;
	const float * ey;
	const float * ez;
	float * w;
	int numPspeeds;
	int * Pspeeds;
	int numMspeeds;
	int * Mspeeds;

	// local partition variables
	int numMySlices,numMyNodes, firstSlice, lastSlice, totalSlices, nnodes;
	float * fEven, *fOdd, * u_bc;
	float rho_out;
	int nd_m, nd_p;
	int * snl, * inl, * onl;

        set<int> myGlobalNodes;
        map<int,int> globalToLocal;

	// MPI communication info
	int firstNdm, lastNdm, firstNdp, lastNdp, numHALO;
	
	// data buffers
	float * ghost_in_m, * ghost_out_m, * ghost_in_p, * ghost_out_p;

	// visualization dump count variable
	int vtk_ts;
	// visualization file offset
	int offset, numEntries;
	float * rho_l, *ux_l, *uy_l, *uz_l;

	TurbineChannel3D(const int rank, const int size);
  ~TurbineChannel3D();
  void write_data(MPI_Comm comm, bool isEven);
  void take_lbm_timestep(bool isEven, MPI_Comm comm);
  
  //pragma acc routine seq
  static inline unsigned getIdx(const unsigned nnodes, const unsigned numSpd, unsigned cellIdx, unsigned spdIdx) 
  { 
     //return  cellIdx*numSpd+spdIdx;
     return spdIdx*nnodes+cellIdx;
  }

 private:
  void read_input_file(const string input_file);
  void initialize_lattice_data();
  void initialize_local_partition_variables();
  void initialize_mpi_buffers();
  void D3Q15_process_slices(bool isEven, const int firstSlice,
		  const int lastSlice);
  void stream_out_collect(bool isEven,const int z_start,float * buff_out, 
      const int numStreamSpeeds, const int * streamSpeeds);
  void stream_in_distribute(bool isEven,const int z_start,
			const float * buff_in, const int numStreamSpeeds,
			const int * streamSpeeds);






};

#endif
