#pragma once
#include<iostream>
#include<opencv2\core\core.hpp>
#include<opencv2\highgui\highgui.hpp>
#include<opencv2\imgproc\imgproc.hpp>
#include<opencv2\imgproc\types_c.h>
#include<vector>
#include"SLIC.hpp"
#include"vibe-background-sequential.hpp"
#include"LBSP\BackgroundSubtractorSuBSENSE.h"
#include"RTSuperpixel.h"
#include<sstream>
#include <fstream>
#include<iomanip>
using namespace std;
using namespace cv;

class SuperPixelBackgroundSubstractor {
public:
	SuperPixelBackgroundSubstractor(String prePath, bool isDown, bool show);
	virtual ~SuperPixelBackgroundSubstractor();
	void generateForegroundSLICVIBE(Mat frame, int frameNumber);
	void generateForegroundSubsense(Mat frame, int frameNumber);
	void generateForegroundGMM(Mat frame, int frameNumber);
	void generateSuperpixel(Mat frame, unsigned short *label, int K, int &realnumber);
	void processWithSLICSuperPixel(Mat frame, Mat foreground, int frameNumber);
	void processWithRTSuperPixel(Mat frame, Mat foreground, int frameNumber, int superpixel);
	vibeModel_Sequential_t *model;
	SLIC slic;
	int  K, spcounta, compactness;
	Mat foreground;
	Mat vibeForeground;
	BackgroundSubtractorSuBSENSE oBGSAlg;
	RTSuperpixel rtSuperpixel;
	Ptr<BackgroundSubtractor> pMOG2;
	String imgnamepre = "bin";
	String imgnametype = ".png";
	String prepath;
	bool isDown;
	bool showFig;
};