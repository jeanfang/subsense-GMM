
#include "SuperPixelBackgroundSubstractor.h"

typedef unsigned int UINT;
using namespace std;


SuperPixelBackgroundSubstractor::SuperPixelBackgroundSubstractor(String prePath, bool down, bool show) {
	model = NULL;
	K = 500;
	spcounta = 10;
	compactness = 20;
	pMOG2 = createBackgroundSubtractorMOG2();
	prepath = prePath;
	isDown = down;
	rtSuperpixel = RTSuperpixel();
	showFig = show;
}

SuperPixelBackgroundSubstractor::~SuperPixelBackgroundSubstractor() {
	model = NULL;
}

void SuperPixelBackgroundSubstractor::generateForegroundSLICVIBE(Mat frame, int frameNumber) {
	Mat gray, suplabel;
	cvtColor(frame, gray, CV_BGR2GRAY);
	if (frameNumber == 1)
	{
		vibeForeground = Mat(gray.rows, gray.cols, CV_8UC1);
		model = (vibeModel_Sequential_t*)libvibeModel_Sequential_New();
		libvibeModel_Sequential_AllocInit_8u_C1R(model, gray.data, gray.cols, gray.rows);
	}
	int width = gray.cols;
	int height = gray.rows;
	int sz = width*height;
	UINT *img = new UINT[sz];

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
			img[i*width + j] = saturate_cast<unsigned  int>(gray.at<uchar>(i, j));
	}
	int *labels = new int[sz];
	int realNumber = 0;

	//slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(img, width, height, labels, realNumber, K, compactness);
	//slic.DrawContoursAroundSegments(img, labels, width, height, 0);

	libvibeModel_Sequential_Segmentation_8u_C1R(model, gray.data, labels, realNumber, vibeForeground.data);
	libvibeModel_Sequential_Update_8u_C1R(model, gray.data, vibeForeground.data);

	/*
	Mat superMask(gray.size(),CV_8UC1, cv::Scalar_<uchar>(255));
	uchar* ptr = superMask.ptr<uchar>(0);
	for (int i = 0; i < height; i++)
	{
	ptr = superMask.ptr<uchar>(i);
	for (int j = 0; j < width; j++) {
	ptr[j] = labels[i*width + j];
	}
	}
	imshow("SLIC", superMask);*/
	if (labels != NULL) delete[] labels;
	if (img != NULL) delete[] img;

	if (isDown) {
		stringstream ss;
		ss << setfill('0') << setw(6) << frameNumber;
		String path = prepath + /*"SlicVibe\\" +*/ imgnamepre + ss.str() + imgnametype;
		imwrite(path, vibeForeground);
	}
	if (showFig)imshow("Vibe-Foreground", vibeForeground);
}

void SuperPixelBackgroundSubstractor::generateForegroundSubsense(Mat frame, int frameNumber) {
	if (frameNumber == 1) {
		foreground.create(frame.size(), CV_8UC1);
		Mat oSequenceROI(frame.size(), CV_8UC1, cv::Scalar_<uchar>(255)); // for optimal results, pass a constrained ROI to the algorithm (ex: for CDnet, use ROI.bmp)
		oBGSAlg.initialize(frame, oSequenceROI);
	}
	oBGSAlg.apply(frame, foreground, double(frameNumber <= 100));
	if (isDown) {
		stringstream ss;
		ss << setfill('0') << setw(6) << frameNumber;
		String path = prepath + /*"Subsense\\" +*/ imgnamepre + ss.str() + imgnametype;
		imwrite(path, foreground);
	}
	if (showFig)imshow("SPBS", foreground);

	//oBGSAlg.process(frame, foreground, double(frameNumber <= 100), labels, realNumber);
	/*processWithSLICSuperPixel(frame, foreground);
	if (isDown) {
	stringstream ss;
	ss << setfill('0') << setw(6) << frameNumber;
	String path = prepath + "SlicSubsense\\" + imgnamepre + ss.str() + imgnametype;
	imwrite(path, foreground);
	}
	imshow("SLIC+Subsense", foreground);*/
	//if (labels != NULL) delete[] labels;
}

void SuperPixelBackgroundSubstractor::generateForegroundGMM(Mat frame, int frameNumber) {
	pMOG2->apply(frame, foreground);
	if (showFig)imshow("GMM-Foreground", foreground);

	if (isDown) {
		stringstream ss;
		ss << setfill('0') << setw(6) << frameNumber;
		String path = prepath + "GMM\\" + imgnamepre + ss.str() + imgnametype;
		imwrite(path, foreground);
	}

	//processWithRTSuperPixel(frame, foreground);
	//imshow("GMM-RT-Foreground", foreground);
}

void SuperPixelBackgroundSubstractor::generateSuperpixel(Mat frame, unsigned short *label, int K, int &realnumber) {
	unsigned char* R, *G, *B;
	int width = frame.cols;
	int height = frame.rows;
	int sz = width*height;

	R = new unsigned char[sz];
	G = new unsigned char[sz];
	B = new unsigned char[sz];
	unsigned int *labnumb = new unsigned int[sz];

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			B[i*width + j] = frame.at<Vec3b>(i, j)[0];
			G[i*width + j] = frame.at<Vec3b>(i, j)[1];
			R[i*width + j] = frame.at<Vec3b>(i, j)[2];
			label[i*width + j] = 0;
		}
	}

	rtSuperpixel.DBscan(R, G, B, height, width, label, K, realnumber, 1, labnumb);

	//DBscan(R, G, B, height, width, label, K, realnumber, 1);
	delete[] R;
	delete[] G;
	delete[] B;
	delete[] labnumb;


	Mat superMask = frame.clone();
	uchar* ptr = superMask.ptr<uchar>(0);
	int superpixelValue = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++) {
			int p = label[i*width + j];
			int count = 0;
			int minx = max(i - 1, 0);
			int maxx = min(i + 1, height - 1);
			int miny = max(j - 1, 0);
			int maxy = min(j + 1, width - 1);
			for (int u = minx; u <= maxx; u++) {
				for (int v = miny; v <= maxy; v++) {
					if (label[u*width + v] != p)
						count++;
					if (count == 2)
						break;
				}
				if (count == 2)
					break;
			}
			if (count == 2) {
				superMask.at<Vec3b>(i, j)[0] = 0;
				superMask.at<Vec3b>(i, j)[1] = 0;
				superMask.at<Vec3b>(i, j)[2] = 0;
			}
		}
	}

	if (showFig)imshow("RTsuperpixel", superMask);
}

void SuperPixelBackgroundSubstractor::processWithSLICSuperPixel(Mat frame, Mat foreground, int frameNumber) {
	Mat gray, suplabel;
	cvtColor(frame, gray, CV_BGR2GRAY);
	int width = gray.cols;
	int height = gray.rows;
	int sz = width*height;
	UINT *img = new UINT[sz];

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
			img[i*width + j] = saturate_cast<unsigned  int>(gray.at<uchar>(i, j));
	}
	int *labels = new int[sz];
	int realNumber = 0;
	slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(img, width, height, labels, realNumber, K, compactness);
	slic.DrawContoursAroundSegments(img, labels, width, height, 0);
	const int *tmplabels = labels;
	//slic.SaveSuperpixelLabels(tmplabels, width, height, "label", "C:\\Users\\YinboYu\\Desktop\\signal statistic\\");
	Mat superMask(gray.size(), CV_8UC1, cv::Scalar_<uchar>(255));
	uchar* ptr = superMask.ptr<uchar>(0);
	for (int i = 0; i < height; i++)
	{
		ptr = superMask.ptr<uchar>(i);
		for (int j = 0; j < width; j++) {
			ptr[j] = img[i*width + j];
		}
	}
	if (showFig)imshow("SLIC", superMask);

	uint8_t *segmentation_map = foreground.data;
	int *numRead = NULL, *isFG = NULL, *numFG = NULL, *numPixel = NULL, *isToken = NULL;
	isFG = (int*)calloc(realNumber, sizeof(int));//is foreground of each superpixel
	numPixel = (int*)calloc(realNumber, sizeof(int));//number of pixel in each superpixel
	numFG = (int*)calloc(realNumber, sizeof(int));//number of foreground pixel in each superpixel
	numRead = (int*)calloc(realNumber, sizeof(int));// the number of read read pixels in each superpixel
	isToken = (int*)calloc(realNumber, sizeof(int));//tag a superpixel 

	memset(numRead, 0, realNumber);
	memset(isFG, 0, realNumber);
	memset(numFG, 0, realNumber);
	memset(numPixel, 0, realNumber);
	memset(isToken, 0, realNumber);

	for (int index = 0; index < width*height; index++)//it is wrong. We cannot directly use the value in each label to tag them. 
	{
		numPixel[labels[index]]++;
	}

	for (int kk = 0; kk<width*height; kk++)
	{
		if (!isToken[labels[kk]]) {
			if (segmentation_map[kk]>0)
				numFG[labels[kk]]++;
			numRead[labels[kk]]++;
			if (double(numRead[labels[kk]]) / double(numPixel[labels[kk]]) > 0.8) {
				isToken[labels[kk]] = 111;
				if (double(numFG[labels[kk]]) / (double)numRead[labels[kk]] > 0.5)
					isFG[labels[kk]] = UCHAR_MAX;
			}
		}
	}

	int kk = 0;
	//imshow("Subsense-Foreground", foreground);
	for (uint8_t *mask = segmentation_map; mask < segmentation_map + (width * height); ++mask, ++kk) {
		/*if (isFG[labels[kk]] == UCHAR_MAX && *mask > 0)
		*mask = UCHAR_MAX;
		else if (isFG[labels[kk]] == UCHAR_MAX)
		*mask = 160;
		else if (*mask > 0)
		*mask = 80;*/
		if (isFG[labels[kk]] == UCHAR_MAX)
			*mask = UCHAR_MAX;
	}

	if (showFig)imshow("Subsense+SLIC", foreground);

	if (isDown) {
		stringstream ss;
		ss << setfill('0') << setw(6) << frameNumber;
		String path = prepath /*+ "SlicSubsense\\" */+ imgnamepre + ss.str() + imgnametype;
		imwrite(path, foreground);
	}

	if (isFG) { free(isFG); isFG = NULL; }
	if (isToken) { free(isToken); isToken = NULL; }
	if (numRead) { free(numRead); numRead = NULL; }
	if (numPixel) { free(numPixel); numPixel = NULL; }
	if (numFG) { free(numFG); numFG = NULL; }


	if (labels != NULL) delete[] labels;
	if (img != NULL) delete[] img;
}

void SuperPixelBackgroundSubstractor::processWithRTSuperPixel(Mat frame, Mat foreground, int frameNumber, int superpixel) {

	unsigned char* R, *G, *B;
	int width = frame.cols;
	int height = frame.rows;
	int sz = width*height;

	R = new unsigned char[sz];
	G = new unsigned char[sz];
	B = new unsigned char[sz];
	unsigned short* labels = new unsigned short[sz];
	unsigned int* labnumb = new unsigned int[sz];
	int realnumber = 0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			B[i*width + j] = frame.at<Vec3b>(i, j)[0];
			G[i*width + j] = frame.at<Vec3b>(i, j)[1];
			R[i*width + j] = frame.at<Vec3b>(i, j)[2];
			labels[i*width + j] = 0;
		}
	}

	rtSuperpixel.DBscan(R, G, B, height, width, labels, superpixel, realnumber, 1, labnumb);

	Mat superMask = frame.clone();
	uchar* ptr = superMask.ptr<uchar>(0);
	int superpixelValue = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++) {
			int p = labels[i*width + j];
			int count = 0;
			int minx = max(i - 1, 0);
			int maxx = min(i + 1, height - 1);
			int miny = max(j - 1, 0);
			int maxy = min(j + 1, width - 1);
			for (int u = minx; u <= maxx; u++) {
				for (int v = miny; v <= maxy; v++) {
					if (labels[u*width + v] != p)
						count++;
					if (count == 2)
						break;
				}
				if (count == 2)
					break;
			}
			if (count == 2) {
				superMask.at<Vec3b>(i, j)[0] = 0;
				superMask.at<Vec3b>(i, j)[1] = 0;
				superMask.at<Vec3b>(i, j)[2] = 0;
			}
		}
	}

	if (showFig)imshow("RTsuperpixel", superMask);


	uint8_t *segmentation_map = foreground.data;

	int *numRead = NULL, *isFG = NULL, *numFG = NULL, *isToken = NULL;
	isFG = (int*)calloc(sz, sizeof(int));//is foreground of each superpixel
	numFG = (int*)calloc(sz, sizeof(int));//number of foreground pixel in each superpixel
	numRead = (int*)calloc(sz, sizeof(int));// the number of read read pixels in each superpixel
	isToken = (int*)calloc(sz, sizeof(int));//tag a superpixel 

	memset(numRead, 0, sz);
	memset(isFG, 0, sz);
	memset(numFG, 0, sz);
	memset(isToken, 0, sz);

	for (int kk = 0; kk<width*height; kk++)
	{
		if (!isToken[labels[kk]]) {
			if (segmentation_map[kk]>0)
				numFG[labels[kk]]++;
			numRead[labels[kk]]++;
			double count = double(numRead[labels[kk]]);
			double N = double(labnumb[labels[kk]]);
			double fore = double(numFG[labels[kk]]);
			if ( count/N > 0.8) {
				isToken[labels[kk]] = 111;
				if ( fore/count > 0.5)
					isFG[labels[kk]] = UCHAR_MAX;
			}
		}
	}

	int kk = 0;
	//imshow("Subsense-Foreground", output);
	for (uint8_t *mask = segmentation_map; mask < segmentation_map + (width * height); ++mask, ++kk) {
		if (*mask > 0) *mask = UCHAR_MAX;
		else if (isFG[labels[kk]] == UCHAR_MAX)
			*mask = 125;
	}
		//if (isFG[labels[kk]] == UCHAR_MAX) *mask = UCHAR_MAX;
	if (showFig)imshow("Subsense+RT", foreground);

	if (isDown) {
		stringstream ss;
		ss << setfill('0') << setw(6) << frameNumber;
		String path = prepath /*+ "RTSubsense\\" */+ imgnamepre + ss.str() + imgnametype;
		imwrite(path, foreground);
	}

	if (numRead) { free(numRead); numRead = NULL; }
	if (isFG) { free(isFG); isFG = NULL; }
	if (numFG) { free(numFG); numFG = NULL; }
	if (isToken) { free(isToken); isToken = NULL; }
	delete[] R;
	delete[] G;
	delete[] B;
	delete[] labels;
	delete[] labnumb;
}