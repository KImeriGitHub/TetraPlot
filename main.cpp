#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<chrono>
#include<numeric>
#include<string>
#include<omp.h> // For parallelization and concurrency
#include<complex>
#include<cmath>

#define cp std::complex<double>

// CALCULATION VARIABLES
const int itermax = 1000;     
const double escTol = 1e30;     
const double approxTol = 1e-30;  

// IMAGE VARIABLES
const char* const name = "tetraPlottest.ppm";
const int XLength = 1920/2;    const int YLength = 1080/2; // Dimension is in pixels
const long double rmi = -5; const long double rma = 3; 

// FORCING EQUISCALED AXIS
/* It workes without equiscaled axis. */
const long double imi = -double(YLength)/XLength*(rma-rmi)/2.0;
const long double ima =  double(YLength)/XLength*(rma-rmi)/2.0;

// SYMMETRY
/* This variable halves the processing time by copying the upper half-sphere.
Only use it if the function below is holomorphic, and imi == -ima . */
const bool YSymmetry = true;

// THREADING PARAMETERS
const int max_threads = omp_get_max_threads();
//const int THREAD_NUM = max_threads-2;
//const int THREAD_NUM = max_threads;
const int THREAD_NUM = 2;

// THE FUNCTION
/* c corresponds to the current point on the complex plane.
   z corresponds to the variable on which we iterate. */
cp f(const cp c, cp z){
	return pow(c,z);
	// return pow(c,-z);
	// return pow(c,-0.5*(z*cos(z)));
	// return pow(c,cos(z));
	// return pow(c,z*sin(z)*0.5);
	// return pow(c,-0.5/z*asin(z)*sin(z));
	// return pow( c, (z*sinh(z)*sin(z)) *-0.5);
	// return cosh(pow(c *(-0.5), z))/c;
	// return pow(c, (z*tan(z))*-0.5);//[-10, 10]
	// return sinh(pow(c*-0.5, z))/c;//[-3.75, 4.25]
	// return pow(c,z*z/c);
	// return pow(c,z*z/c/c);
	// return exp(pow(c,z));
}

// PRE: o is a empty file opened through ifstream in binary
// POST: A ppm file with at most 3*XLength*YLength of 8-bit data-values.
void writeppm(std::ostream& o);

// POST: returns an integer corresponding to the behaviour of cp
int tetracheck(const cp c);

// POST: 3 8-bit values are written on color, corresponding to k.
void coloring(const int k, unsigned char color[3]);

int main(){
	// Measuring elapsed time
	std::chrono::time_point<std::chrono::system_clock> tic, toc;
	tic = std::chrono::system_clock::now();

	// Open file and process
	std::ofstream ppmfile(name, std::ios::binary);
    writeppm(ppmfile);
    ppmfile.close();

	// We copy the upper half-space onto the lower one.
	if (YSymmetry){
		// Open processed file to read.
		std::ifstream ppmInFile(name, std::ios::binary);
		std::string line;

		// Skip first lines
		std::getline(ppmInFile, line); // P6
		std::getline(ppmInFile, line); // Xlength Ylength
		std::getline(ppmInFile, line); // 255

		// Copy the whole picture onto the string line
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// !! THIS MAY DEMAND TOO MUCH SPACE IN THE RAM !!
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		std::getline(ppmInFile, line); 

		// The following creates a pointer to the string
		const char* linePtr = line.c_str();
		ppmInFile.close();

		// Open processed file to write.
		std::ofstream ppmfile(name, std::ios::binary | std::ios::app);
		// We need to copy the image upside down.
		const char* linePartPtr; 
		for(int i = 0; i<YLength/2;i++){
			// We set linePartPtr at the corresponding row-index and copy that row.
			linePartPtr = linePtr+(YLength/2-i-1)*3*XLength;
			ppmfile.write(linePartPtr, 3*XLength);
		}
		ppmfile.close();
	}

	// Final terminal output
	toc = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsedTime = toc - tic;
	int elaHours = int(elapsedTime.count()/3600);
	int elaMinutes = int((elapsedTime.count() - elaHours*3600)/60);
	std::cout<<"\n Program ended in " << elaHours   << " hours "
									  << elaMinutes << " minutes "
									  << elapsedTime.count() - elaHours*3600 - elaMinutes*60 << " seconds.";

	return 0;
}

void writeppm(std::ostream& o){
	// Initial information on a ppm file
    o << "P6" << std::endl
      << XLength << " " << YLength << std::endl
      << "255" << std::endl;

	// We save the elapsed time for every row in timevec
	std::vector<double> timeVec;
	timeVec.push_back(0);
    std::chrono::time_point<std::chrono::system_clock> timeStart, timeEnd;

	// YLengthMod is for symmetric cases
	int YLengthMod = YLength;
	if(YSymmetry) YLengthMod/=2;

	/////////////////////////////////////
	// Loop over every row			   //
	/////////////////////////////////////
	int progmarker=0; //for marking progress; it is an integer between 0 and 100.
	for(int j=0; j < YLengthMod; ++j){
		unsigned char color_arr[3*XLength];  // coloring array

		// im: Current imaginary value of the complex number
        const long double im = ima - (ima -imi) *j /(YLength -1);

		// Informations on the terminal: 
		double timeVecSum = std::reduce(timeVec.begin(), timeVec.end());
		int timeVecSize = timeVec.size();
		if(double(progmarker)<=(100*j/(YLengthMod-1))){
			double estSeconds = (timeVecSum / std::max((double)progmarker,double(1)) * (100-(double)progmarker));
			int estHours = int(estSeconds/3600);
			int estMinutes = int((estSeconds - estHours*3600)/60);
			std::cout<< progmarker <<" %   "<<
				 "Estimated time left: " << estHours << " Hours " 
				  				    	 << estMinutes << " Minutes "
								    	 << estSeconds - estHours*3600 - estMinutes*60 << " Seconds " <<std::endl;
			++progmarker;
		}

		///////////////////////////////////
		// Main Loop					 //
		///////////////////////////////////
		timeStart = std::chrono::system_clock::now();

		// Threading prerequisites
		// For more information we refer to 
		// https://www.openmp.org/resources/tutorials-articles/
		omp_set_num_threads(THREAD_NUM); // set number of threads in "parallel" blocks
		int nthreads;

		// DO NOT INDENT
		#pragma omp parallel 
		{
		int i, id,nthrds;
		id = omp_get_thread_num();
		nthrds=omp_get_num_threads();
		if(id==0) nthreads = nthrds; // Main thread
		for (i=id; i < XLength; i+=nthrds) // Every thread gets its own for loop.
		{
			// re: Current real value of the complex number 
			const long double re = rmi + (rma - rmi) * i / (XLength - 1);
			const cp c(re, im);

			// Assign the 3 8-bit values to color_loc
			unsigned char color_loc[3]; 
			coloring(tetracheck(c), color_loc);

			//Assign local coloring to global color array.
			color_arr[i*3+0] = color_loc[0];
			color_arr[i*3+1] = color_loc[1];
			color_arr[i*3+2] = color_loc[2];
		}
		} // End Parallelisation

		// Writing whole row at once onto the file
		o.write((const char*)color_arr,3*XLength);

		// Measuring elapsed time
		timeEnd = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = timeEnd - timeStart;
		timeVec.push_back(elapsed_seconds.count());
    }
}

int tetracheck(const cp c){
	// Given c, we apply f(c, . ) iteratively,
	// that is f(c, f(c, ... c, f(c, f(c, 0))...)), 
	// until we either see a divergence, a convergence or a cycle.
	// On divergence, it returns 0.
	// On convergence it returns 1.
	// On cycle it returns the cycle length.
	
	cp z(0,0);			//Initial value of the iteration

	// zVec stores each complex value of the iteration
	std::vector<cp> zVec(itermax, z);
	
    for(int i=0;i<itermax;i++){ // THIS FOR-LOOP HAS TWO ESCAPE LINES

		z=f(c,z);       // Here is the function!
		zVec[i] = z;

		// Check whether it escapes
		if(abs(z)>escTol) {
			return i;
		} 

		// If it did not escape check whether it cycles
		for (int j=0;j<i;j++){
			if (abs(zVec[j]-z) < approxTol)
				return i-j;
		}
	}
	return 0;
}

void coloring(const int k, unsigned char color[]){
	if (k==0)      {color[0] = 0;   color[1] = 0;   color[2] = 0;}
	else if (k==1) {color[0] = 255; color[1] = 0;   color[2] = 0;}  //red
	else if (k==2) {color[0] = 120; color[1] = 15;  color[2] = 15;} //darkred
	else if (k==3) {color[0] = 0;   color[1] = 0;   color[2] = 215;}//blue
	else if (k==4) {color[0] = 213; color[1] = 105; color[2] = 0;}  //orange
	else if (k==5) {color[0] = 75;  color[1] = 4;   color[2] = 162;}//magenta
	else if (k==6) {color[0] = 0;   color[1] = 255; color[2] = 255;}//cyan
	else if (k==7) {color[0] = 255; color[1] = 255; color[2] = 255;}//white
	else{color[0] = int(255. / log(double(k - 5)));
	     color[1] = int(255. / log(double(k - 5)));
	     color[2] = int(255. / log(double(k - 5)));
	}
}