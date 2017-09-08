#include <iostream>
#include <string>
#include <valarray>
#include <vector>
#include <iterator>
#include <thread>
#include <CCfits/CCfits>


using namespace std;
//using namespace CCfits;


/*
 ********* Parameters for threads
 */
struct ParamsThread {
	unsigned int index;		// Index of thread
	unsigned int startY;	// Global index, from where to start computations and where to write the result
	unsigned int endY;		// Size of the computations for thread
};


class Averager {

private:
	string fileName;
	int startX,
		endX,
		startY,
		endY,
		startZ,
		endZ,
		numThreads = 2;
	double cdelt;
	valarray<float> data, res;
	vector<unsigned int> dataDimentions;

public:
	Averager(char* fileName) {
		readImage(fileName);
		this->res.resize(this->dataDimentions[0] * this->dataDimentions[1]);
	}

	Averager(char* fileName, int startX, int endX, int startY, int endY, int startZ, int endZ, double cdelt) {
		this->fileName = fileName;
		this->startX = startX;
		this->endX = endX;
		this->startY = startY;
		this->endY = endY;
		this->startZ = startZ;
		this->endZ = endZ;
		this->cdelt = cdelt;

		cout << fileName << endl;
		readImage(fileName);
		this->res.resize(this->dataDimentions[0] * this->dataDimentions[1]);
	}

	Averager(char* fileName, int startX, int endX, int startY, int endY, int startZ, int endZ, double cdelt, int numThreads) {
		this->fileName = fileName;
		this->startX = startX;
		this->endX = endX;
		this->startY = startY;
		this->endY = endY;
		this->startZ = startZ;
		this->endZ = endZ;
		this->cdelt = cdelt;
		this->numThreads = numThreads;

		cout << this->numThreads << endl;
		readImage(fileName);
		this->res.resize(this->dataDimentions[0] * this->dataDimentions[1]);
	}

	Averager(int startX, int endX, int startY, int endY, int startZ, int endZ, double cdelt) {
			this->startX = startX;
			this->endX = endX;
			this->startY = startY;
			this->endY = endY;
			this->startZ = startZ;
			this->endZ = endZ;
			this->cdelt = cdelt;

			cout << cdelt << endl;
		}

	Averager() {}

	void setVariables(int startX, int endX, int startY, int endY, int startZ, int endZ, double cdelt, int numThreads) {
		this->startX = startX;
		this->endX = endX;
		this->startY = startY;
		this->endY = endY;
		this->startZ = startZ;
		this->endZ = endZ;
		this->cdelt = cdelt;
		this->numThreads = numThreads;
	}

	void readImage(char* fileName) {

	/*	if( info )
		  FITS::setVerboseMode(true);*/

	  try {
		   auto_ptr<CCfits::FITS> pInfile( new CCfits::FITS( fileName, CCfits::Read, true ) );
		   CCfits::PHDU& image = pInfile->pHDU();

		   // read all user-specifed, coordinate, and checksum keys in the image
		  // image.readAllKeys();
		   image.read( this->data );

		   int size =  image.axes();
		   this->dataDimentions.resize( size );

		   for(int i = 0; i < size; ++i )
			   this->dataDimentions[i] = image.axis( i );

		} catch (CCfits::FitsException&) {
		  cerr << " Fits Exception Thrown by FitsImage class \n";
		  cerr << " Fits file name : " << fileName << endl;
		}
		catch (std::bad_alloc& ba) {
		   cerr << "bad_alloc caught: " << ba.what() << endl;
		   cerr << "Not enough space in RAM to load FITS file" << endl;
		   exit(1);
		 }
	}

	vector<unsigned int> getDataDimentions() {

		return this->dataDimentions;
	}



	void runThreads(unsigned int index, unsigned int startY, unsigned int endY) {

		cout << "Thread #" << index << " will process: " << startY << ", " << endY << endl;

		unsigned int step = this->dataDimentions[0] * this->dataDimentions[1];
		unsigned int count = this->endZ - this->startZ;
		unsigned int offset = this->startZ * step;

		for(unsigned int i = startY; i < endY; ++i) {
			unsigned int realInd = i * this->dataDimentions[0];

			for(unsigned int j = 0; j < this->dataDimentions[0] ; ++j) {

				unsigned int resIj = realInd + j;
				unsigned int start = offset + resIj;

				valarray<float> tmp = this->data[slice(start, count, step)];
				this->res[resIj] = tmp.sum() * this->cdelt;
			}
		}

	}


	vector<float> countAverage() {


		cout << "data size: " << this->data.size() << " result size: " << this->res.size() << endl;

		vector<thread> threads( this->numThreads );
		vector<ParamsThread> paramsThreads( this->numThreads );

		unsigned int jobSize = this->dataDimentions[1] / this->numThreads;
		unsigned int remainder = this->dataDimentions[1] % this->numThreads;

		paramsThreads[0].index = 0;
		paramsThreads[0].startY = 0;
		paramsThreads[0].endY = paramsThreads[0].startY + (0 < remainder ? jobSize + 1 : jobSize);

		for (unsigned int i = 1; i < this->numThreads; ++i) {
			paramsThreads[i].index = i;
			paramsThreads[i].startY = paramsThreads[i - 1].endY; // + paramsThreads[i - 1].startY;
			paramsThreads[i].endY = paramsThreads[i].startY + (i < remainder ? jobSize + 1 : jobSize);
		}

		for (unsigned int i = 0; i < remainder; ++i) {
			++paramsThreads[i].endY;
			if( i < ( paramsThreads.size() - 1 ) )
				++paramsThreads[i + 1].startY;
		}

		for (unsigned int i = 0; i < this->numThreads; ++i) {
			threads[i] = thread(&Averager::runThreads, this, paramsThreads[i].index, paramsThreads[i].startY, paramsThreads[i].endY);
		}

		for (unsigned int i = 0; i < this->numThreads; ++i) {
			threads[i].join();
		}

		vector<float> v(begin(this->res), end(this->res));

		return v;
	}

};



