#include <iostream>
#include <opencv2\core\core.hpp>
#include <opencv2\video\background_segm.hpp>
#include <opencv2\videoio.hpp>
#include <opencv2\imgcodecs.hpp>
#include <opencv2\video.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\imgproc\types_c.h>
#include <vector>
#include "SLIC.hpp"
#include "vibe-background-sequential.hpp"
#include "SuperPixelBackgroundSubstractor.h"
#include "LBSP\BackgroundSubtractorSuBSENSE.h"
#include <sstream>
#include <iomanip>
typedef unsigned int UINT;
using namespace std;
using namespace cv;


int main()
{
	Mat frame;
	VideoCapture capture;
	//String vidoePrePath = "D:\\BJJM712\\BJJM_dataset_all\\dynamicBackground\\";
	String vidoePrePath = "..\\dataset\\";
	String videos[6] = { "fall", "streetLight", "office", "intermittentPan", "highway", "canoe" };
	
	String video = "\\input\\in%06d.jpg";
	String results = "_results\\";

	String outPath = "..\\Data\\BGS";
	for (int i = 0; i < 1; i++) {
		//VideoCapture sequence(vidoePrePath + videos[i] + video);
		capture.open(vidoePrePath + videos[i] + video);
		if (capture.isOpened()) {
			cout << "Open and process the video of " << videos[i] << endl;
			String ss = "\\"+videos[i];
			SuperPixelBackgroundSubstractor spbgs = SuperPixelBackgroundSubstractor(outPath + ss + results, true, true);
			RTSuperpixel rt = RTSuperpixel();
			for (int k = 0;; ++k) {
				if (!capture.read(frame)) {
					break;
				}
				imshow("input", frame);
				int pixels = frame.cols*frame.rows;
				waitKey(1);
				cout << "Process the " << k<<" frame" << endl;
				spbgs.generateForegroundSubsense(frame, k + 1);

				//GMM
				spbgs.generateForegroundGMM(frame, k + 1);
				//spbgs.generateForegroundSLICVIBE(frame, k + 1);

				if (k == 325)
					cout<<"";
				
			}
		}
	}

	
	return 0;
}





