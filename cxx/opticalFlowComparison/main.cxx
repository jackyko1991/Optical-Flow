#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkOpenCVImageBridge.h"
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkOpenCVImageBridge.h>
#include <itkCastImageFilter.h>
#include <itkMinimumProjectionImageFilter.h>
#include <itkSubtractImageFilter.h>

#include <OpticalFlow.h>
#include <ImageIO.h>
#include <direct.h>

#include <ctime>

#include <QFileInfo>
#include <QString>
#include <QTextStream>

// includes from OpenCV
#include "cv.h"

#if CV_VERSION_MAJOR > 2
#include "opencv2/opencv.hpp" // cv::imwrite
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/legacy/legacy.hpp>
#endif

using namespace std;
using namespace cv;

typedef itk::Image<float, 3> Image3DType;
typedef itk::Image<unsigned char, 2> Image2DUCType;
typedef itk::Image<float, 2> Image2DFloatType;
typedef itk::ImageFileReader< Image3DType >  ReaderType;
typedef itk::ImageRegionIterator<Image3DType> Iterator3DType;
typedef itk::ImageRegionIterator<Image2DFloatType> Iterator2DType;
typedef itk::ImageFileWriter<Image2DUCType> WriterType;
typedef itk::RescaleIntensityImageFilter<Image2DFloatType, Image2DUCType> RescaleFilterType;
typedef itk::CastImageFilter< Image2DFloatType, Image2DUCType > CastImageFilterType;
typedef itk::MinimumProjectionImageFilter<Image3DType, Image2DFloatType> minIPFilterType;
typedef itk::SubtractImageFilter<Image2DFloatType> SubtractImageFilterType;

void ExtractImageFrom3D(Image3DType::Pointer Image3D, Image2DFloatType::Pointer Image2D, unsigned int slice)
{
	// 3D region info
	Image3DType::RegionType region3D;
	Image3DType::IndexType index3D;
	index3D[0] = 0;
	index3D[1] = 0;
	index3D[2] = slice;
	Image3DType::SizeType size3D;
	size3D[0] = Image3D->GetLargestPossibleRegion().GetSize()[0];
	size3D[1] = Image3D->GetLargestPossibleRegion().GetSize()[1];
	size3D[2] = 1;
	region3D.SetIndex(index3D);
	region3D.SetSize(size3D);

	// 2D region info
	Image2DFloatType::RegionType region2D;
	Image2DFloatType::IndexType start;
	start[0] = 0;
	start[1] = 0;

	Image2DFloatType::SizeType size;
	size[0] = Image3D->GetLargestPossibleRegion().GetSize()[0];
	size[1] = Image3D->GetLargestPossibleRegion().GetSize()[1];

	region2D.SetSize(size);;
	region2D.SetIndex(start);

	// iterate over the slice
	Iterator3DType it3D(Image3D, region3D);
	Iterator2DType it2D(Image2D, region2D);
	while (!it3D.IsAtEnd())
	{
		it2D.Set(it3D.Value());
		++it3D;
		++it2D;
	}
}

void ConvertToCVImage(Image2DFloatType::Pointer ITK2DImageInput, cv::Mat& cvImageOutput)
{
	// cast and rescale image from float to unsigned char
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(ITK2DImageInput);
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();

	// cast image from itk to opencv
	cv::Mat img = itk::OpenCVImageBridge::ITKImageToCVMat< Image2DUCType >(rescaleFilter->GetOutput());

	cvImageOutput = img.clone();

	cv::Point2f ptCp(img.cols*0.5, img.rows*0.5);
	cv::Mat M = cv::getRotationMatrix2D(ptCp, 180, 1.0);
	cv::warpAffine(img, cvImageOutput, M, img.size(), cv::INTER_CUBIC); //Nearest is too rough, 
}

void Extract2ConsecutiveCVImages(Image3DType::Pointer itkImage, cv::Mat& slice0, cv::Mat& slice1, int slice0Num)
{
	Image2DFloatType::RegionType region2D;
	Image2DFloatType::IndexType start;
	start[0] = 0;
	start[1] = 0;

	Image2DFloatType::SizeType size;
	size[0] = itkImage->GetLargestPossibleRegion().GetSize()[0];
	size[1] = itkImage->GetLargestPossibleRegion().GetSize()[1];

	region2D.SetSize(size);;
	region2D.SetIndex(start);

	Image2DFloatType::Pointer frame0 = Image2DFloatType::New();
	frame0->SetRegions(region2D);
	frame0->Allocate();

	Image2DFloatType::Pointer frame1 = Image2DFloatType::New();
	frame1->SetRegions(region2D);
	frame1->Allocate();

	ExtractImageFrom3D(itkImage, frame0, slice0Num);
	ExtractImageFrom3D(itkImage, frame1, slice0Num+1);

	// convert itk image to opencv image
	ConvertToCVImage(frame0, slice0);
	ConvertToCVImage(frame1, slice1);
}

double OpticalFlowFarneback(cv::Mat& frame0, cv::Mat& frame1, cv::Mat& flow, double pyr_scale, int levels, int winSize, int iter_n, int poly_n, double poly_sigma, int flags, bool displayTiming)
{
	cv::Mat vel;

	clock_t start_t = clock();
	cv::calcOpticalFlowFarneback(frame0, frame1, vel, pyr_scale, levels, winSize, iter_n, poly_n, poly_sigma, flags);
	clock_t end_t = clock();
	double time = (double)(end_t - start_t) / CLOCKS_PER_SEC * 1.0;
	if (displayTiming)
		cout << "Computation time for Farneback method = " << time << " second(s)." << endl;
	
	flow.create(vel.rows, vel.cols, CV_8UC1);
	for (int i = 0; i < vel.rows; ++i)
	{
		for (int j = 0; j < vel.cols; ++j)
		{
			Vec2f flow_at_point = vel.at<Vec2f>(i, j);
			float fx = flow_at_point[0];
			float fy = flow_at_point[1];

			flow.at<unsigned char>(i, j) = sqrt(fx*fx + fy*fy)*50; // scaling is to fit jpg range
		}
	}
	return time;
}

double OpticalFlowLiu(cv::Mat& frame0, cv::Mat& frame1, cv::Mat& flow, double alpha, double ratio, int minWidth, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations, bool displayTiming)
{
	DImage vx, vy, warpI2, Im0, Im1;
	Im0.loadCVImage(frame0);
	Im1.loadCVImage(frame1);

	clock_t start_t = clock();
	OpticalFlow::Coarse2FineFlow(vx, vy, warpI2, Im0, Im1, alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations, nSORIterations);
	clock_t end_t = clock();
	double time = (double)(end_t - start_t) / CLOCKS_PER_SEC * 1.0;
	if (displayTiming)
		cout << "Computation time for Liu method = " << time << " second(s)." << endl;

	// get the velocity magnitude image
	DImage mag, vx2, vy2;
	vx2.copy(vx);
	vx2.square();
	vy2.copy(vy);
	vy2.square();
	mag.Add(vx2, vy2);
	mag.Sqrt();

	DImage scaleMatrix;
	scaleMatrix.copy(mag);
	scaleMatrix.setValue(50.0/255.0); // scaling is to fit jpg range
	mag.Multiplywith(scaleMatrix);

	// return the image in CV type
	mag.getCVImage(flow, ImageIO::standard);

	return time;
}

double OpticalFlowHS(cv::Mat& frame0, cv::Mat& frame1, cv::Mat& flow, int max_ter, double epsilon, double lambda, bool displayTiming)
{
	CvTermCriteria IterCriteria;
	IterCriteria.type = CV_TERMCRIT_ITER | CV_TERMCRIT_EPS;
	IterCriteria.max_iter = max_ter;
	IterCriteria.epsilon = epsilon;

	cv::Mat velx = cv::Mat::zeros(frame0.size().width, frame0.size().height, CV_32FC1);
	cv::Mat vely = cv::Mat::zeros(frame0.size().width, frame0.size().height, CV_32FC1);

	// HS method is old optical flow method in openCV, need to use Ipl Image input
	IplImage frame0arr = frame0;
	IplImage frame1arr = frame1;
	IplImage velxArr = velx;
	IplImage velyArr = vely;

	/* Do optical flow computation */
	clock_t start_t = clock();
	cvCalcOpticalFlowHS(&frame0arr, &frame1arr, 0, &velxArr, &velyArr, lambda, IterCriteria);
	clock_t end_t = clock();
	double time = (double)(end_t - start_t) / CLOCKS_PER_SEC * 1.0;
	if (displayTiming)
		cout << "Computation time for Horn Schunck method = " << time << " second(s)." << endl;

	flow.create(velx.rows, velx.cols, CV_8UC1);
	for (int i = 0; i < velx.rows; ++i)
	{
		for (int j = 0; j < velx.cols; ++j)
		{
			float fx = CV_IMAGE_ELEM(&velxArr, float, i, j);
			float fy = CV_IMAGE_ELEM(&velyArr, float, i, j);

			flow.at<unsigned char>(i, j) = sqrt(fx*fx + fy*fy) * 50; // scaling is to fit jpg range
		}
	}

	return time;
}

double OpticalFlowLK(cv::Mat& frame0, cv::Mat& frame1, cv::Mat& flow, int winSizeX, int winSizeY, bool displayTiming)
{
	cv::Mat velx = cv::Mat::zeros(frame0.size().width, frame0.size().height, CV_32FC1);
	cv::Mat vely = cv::Mat::zeros(frame0.size().width, frame0.size().height, CV_32FC1);

	// LK method is old optical flow method in openCV, need to use Ipl Image input
	IplImage frame0arr = frame0;
	IplImage frame1arr = frame1;
	IplImage velxArr = velx;
	IplImage velyArr = vely;

	CvSize size;
	size.height = winSizeX;
	size.width = winSizeY;

	/* Do optical flow computation */
	clock_t start_t = clock();
	cvCalcOpticalFlowLK(&frame0arr, &frame1arr, size, &velxArr, &velyArr);
	clock_t end_t = clock();
	double time = (double)(end_t - start_t) / CLOCKS_PER_SEC * 1.0;
	if (displayTiming)
		cout << "Computation time for Horn Schunck method = " << time << " second(s)." << endl;

	flow.create(velx.rows, velx.cols, CV_8UC1);
	for (int i = 0; i < velx.rows; ++i)
	{
		for (int j = 0; j < velx.cols; ++j)
		{
			float fx = CV_IMAGE_ELEM(&velxArr, float, i, j);
			float fy = CV_IMAGE_ELEM(&velyArr, float, i, j);

			flow.at<unsigned char>(i, j) = sqrt(fx*fx + fy*fy) * 50; // scaling is to fit jpg range
		}
	}

	return time;
}

void VesselSegmentationByMinIP(Image3DType::Pointer image, Image2DFloatType::Pointer seg)
{
	minIPFilterType::Pointer minIPFilter = minIPFilterType::New();
	minIPFilter->SetInput(image);
	minIPFilter->Update();

	//Extract 1st slice for subtraction
	Image2DFloatType::RegionType region2D;
	Image2DFloatType::IndexType start;
	start[0] = 0;
	start[1] = 0;

	Image2DFloatType::SizeType size;
	size[0] = image->GetLargestPossibleRegion().GetSize()[0];
	size[1] = image->GetLargestPossibleRegion().GetSize()[1];

	region2D.SetSize(size);;
	region2D.SetIndex(start);

	Image2DFloatType::Pointer slice0 = Image2DFloatType::New();
	slice0->SetRegions(region2D);
	slice0->Allocate();

	ExtractImageFrom3D(image, slice0, 0);
	SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New();
	subtractFilter->SetInput1(slice0);
	subtractFilter->SetInput2(minIPFilter->GetOutput());
	subtractFilter->Update();

	seg->Graft(subtractFilter->GetOutput());
}

void ExtractLargestConnectedComponent(const cv::Mat &src, cv::Mat& binaryLCC)
{
	// binarize the image
	cv::threshold(src, binaryLCC, 0, 1, 0);

	// Fill the label_image with the blobs
	// 0  - background
	// 1  - unlabelled foreground
	// 2+ - labelled foreground

	binaryLCC.convertTo(binaryLCC, CV_32SC1);
	cv::Mat label_image = binaryLCC.clone();

	int label_count = 2; // starts at 2 because 0,1 are used already
	int max_area = 0;

	int LCClabelIdx = 1;
	for (int y = 0; y < label_image.rows; y++) {
		int *row = (int*)label_image.ptr(y);
		for (int x = 0; x < label_image.cols; x++) {
			if (row[x] != 1)
			{
				continue;
			}
			cv::Rect rect;
			cv::floodFill(label_image, cv::Point(x, y), label_count, &rect, 0, 0, 4);

			// count the size of area in each label
			int area = 0;

			for (int i = rect.y; i < (rect.y + rect.height); i++) 
			{
				int *row2 = (int*)label_image.ptr(i);
				for (int j = rect.x; j < (rect.x + rect.width); j++) 
				{
					if (row2[j] != label_count) 
					{
						continue;
					}
					area++;
				}
			}

			if (area > max_area)
			{
				LCClabelIdx = label_count;
				max_area = area;
			}

			label_count++;
		}
	}

	binaryLCC = cv::Mat::zeros(binaryLCC.size(), CV_8U);
	// image thresholding to extract LCC
	for (int y = 0; y < label_image.rows; y++)
	{
		int *row = (int*)label_image.ptr(y);
		for (int x = 0; x < label_image.cols; x++)
		{
			if (row[x] == LCClabelIdx)
			{
				binaryLCC.at<unsigned char>(y, x) = 255;
			}
		}
	}
	binaryLCC.convertTo(binaryLCC, 0); // the image type is modified after LCC extraction
	cv::threshold(binaryLCC, binaryLCC, 0, 255, 0);
}

int main(int argc, char* argv[])
{
	if (argc < 3 || argc >5)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0];
		std::cerr << " <InputFileName> <OutputFolderName> <StartSlice> <EndSlice>, if slice number are not specified, all slices will be computed";
		std::cerr << std::endl;
		return EXIT_FAILURE;
	}

	QFileInfo check_input_file(argv[1]);
	if (!(check_input_file.exists() && check_input_file.isFile()))
	{
		cerr << "Input file not exist" << endl;
		return EXIT_FAILURE;
	}

	QFileInfo check_output_folder(argv[2]);
	if (!check_output_folder.exists())
	{
		mkdir(argv[2]);
	}
	string saveFolder(argv[2]);

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);
	cout << "Processing data " << argv[1] << endl;

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & error)
	{
		std::cerr << "Error: " << error << std::endl;
		return EXIT_FAILURE;
	}

	// loop over selected slices, 1st slice marked as 1
	int start = 1;
	int end = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1;

	if (argc == 5)
	{
		start = strtol(argv[3], NULL, 10);
		end = strtol(argv[4], NULL, 10);
	}
	else
	{
		cout << "Start and end slice not setted, all slices will be computed." << endl;
	}

	if (start > end)
	{
		cerr << "Start slice number should be smaller than end slice number" << endl;
		return EXIT_FAILURE;
	}

	if (start < 1 || end > reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1)
	{
		cerr << "Selected slices exceed image range. Please note that the maximum end slice should be one slice less than total number of slices. Total number of slice is " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] << endl;
		return EXIT_FAILURE;
	}

	cout << "Start slice: " << start << endl;
	cout << "End slice: " << end << endl;

	// vessel segmentation by subtracting 1st slice from minIP image
	Image2DFloatType::Pointer vessel = Image2DFloatType::New();
	cv::Mat vesselCV, vesselMaskCV, vesselMaskCV_;

	VesselSegmentationByMinIP(reader->GetOutput(), vessel);
	ConvertToCVImage(vessel, vesselCV);
	/* 0: Binary
	1: Binary Inverted
	2: Threshold Truncated
	3: Threshold to Zero
	4: Threshold to Zero Inverted
	*/
	cv::threshold(vesselCV, vesselMaskCV, 45, 256, 0);
	int morph_elem = 2; // this declare for a ellipse shpae kernel
	int morph_size = 2; 
	Mat element = getStructuringElement(morph_elem, Size(2 * morph_size + 1, 2 * morph_size + 1), Point(morph_size, morph_size)); // kernel for image closing
	cv:morphologyEx(vesselMaskCV, vesselMaskCV, 3, element); // morphological close for noise reduction
	ExtractLargestConnectedComponent(vesselMaskCV, vesselMaskCV);
	cv::bitwise_not(vesselMaskCV, vesselMaskCV_);

	double timeFarneback = 0;
	double timeLiu = 0;
	double timeHS = 0;
	double timeLK = 0;
	int count = 0;

	for (int i = start-1; i <end; i++)
	{
		cout << "Computing slice " << (i + 1) << endl;
		//Extract 2 consecuctive images for 2 frame optical flow calculation
		cv::Mat frame0CV;
		cv::Mat frame1CV;
		Extract2ConsecutiveCVImages(reader->GetOutput(), frame0CV, frame1CV,i);

		// optical flow calculation
		cv::Mat flowFraneback, flowLiu, flowHS, flowLK;
		timeFarneback = timeFarneback + OpticalFlowFarneback(frame0CV,frame1CV,flowFraneback,0.75,3,15,20,5,1.2,0,false);
		timeLiu = timeLiu + OpticalFlowLiu(frame0CV, frame1CV, flowLiu, 0.0006, 0.55, 256, 1, 1 , 20, false);
		timeHS = timeHS + OpticalFlowHS(frame0CV, frame1CV, flowHS, 500, 0.0001, 0.012, false);
		timeLK = timeLK + OpticalFlowLK(frame0CV, frame1CV, flowLK, 15, 15, false);

		// color mapping
		cv::Mat flowFranebackC, flowLiuC, flowHSC, flowLKC;
		cv::applyColorMap(flowFraneback, flowFranebackC, cv::COLORMAP_JET);
		cv::applyColorMap(flowLiu, flowLiuC, cv::COLORMAP_JET);
		cv::applyColorMap(flowHS, flowHSC, cv::COLORMAP_JET);
		cv::applyColorMap(flowLK, flowLKC, cv::COLORMAP_JET);

		// image masking with vessel mask
		cv::Mat flowFranebackCMask, flowLiuCMask, flowHSCMask, flowLKCMask;
		flowFranebackC.copyTo(flowFranebackCMask, vesselMaskCV);
		flowLiuC.copyTo(flowLiuCMask, vesselMaskCV);
		flowHSC.copyTo(flowHSCMask, vesselMaskCV);
		flowLKC.copyTo(flowLKCMask, vesselMaskCV);

		// add the flow to original image
		cv::Mat frame0CVMask;
		frame0CV.copyTo(frame0CVMask, vesselMaskCV_);

		// image need to first convert to RGB space
		cv::cvtColor(frame0CVMask, frame0CVMask, CV_GRAY2RGB);

		double alpha = 0.5; //blending weight
		cv::Mat flowFranebackOverlay, flowLiuOverlay, flowHSOverlay, flowLKOverlay;
		cv::add(frame0CVMask, flowFranebackCMask, flowFranebackOverlay);
		cv::add(frame0CVMask, flowLiuCMask, flowLiuOverlay);
		cv::add(frame0CVMask, flowHSCMask, flowHSOverlay);
		cv::add(frame0CVMask, flowLKCMask, flowLKOverlay);

		// save the images
		string OriginalSavePath;
		string FarnebackSavePath;
		string LiuSavePath;
		string HSSavePath;
		string LKSavePath;

		if (i < 9)
		{
			OriginalSavePath = saveFolder + "\\Original_00" + to_string(i + 1) + ".jpg";
			FarnebackSavePath = saveFolder + "\\Farneback_00" + to_string(i + 1) +".jpg";
			LiuSavePath = saveFolder + "\\Liu_00" + to_string(i + 1) +".jpg";
			HSSavePath = saveFolder + "\\HornSchunck_00" + to_string(i + 1) +".jpg";
			LKSavePath = saveFolder + "\\LucasKanade_00" + to_string(i + 1) +".jpg";
		}
		else if (i < 99)
		{
			OriginalSavePath = saveFolder + "\\Original_0" + to_string(i + 1) + ".jpg";
			FarnebackSavePath = saveFolder + "\\Farneback_0" + to_string(i + 1) + ".jpg";
			LiuSavePath = saveFolder + "\\Liu_0" + to_string(i + 1) + ".jpg";
			HSSavePath = saveFolder + "\\HornSchunck_0" + to_string(i + 1) + ".jpg";
			LKSavePath = saveFolder + "\\LucasKanade_0" + to_string(i + 1) + ".jpg";
		}
		else
		{
			OriginalSavePath = saveFolder + "\\Original_" + to_string(i + 1) + ".jpg";
			FarnebackSavePath = saveFolder + "\\Farneback_" + to_string(i + 1) + ".jpg";
			LiuSavePath = saveFolder + "\\Liu_" + to_string(i + 1) + ".jpg";
			HSSavePath = saveFolder + "\\HornSchunck_" + to_string(i + 1) + ".jpg";
			LKSavePath = saveFolder + "\\LucasKanade_" + to_string(i + 1) + ".jpg";
		}
		
		cv::imwrite(OriginalSavePath, frame0CV);
		cv::imwrite(FarnebackSavePath, flowFranebackOverlay);
		cv::imwrite(LiuSavePath, flowLiuOverlay);
		cv::imwrite(HSSavePath, flowHSOverlay);
		cv::imwrite(LKSavePath, flowLKOverlay);
		count++;
	}

	cout << "Average time for Farneback method = " << timeFarneback/count <<"second(s)" << endl;
	cout << "Average time for Liu method = " << timeLiu / count << "second(s)" << endl;
	cout << "Average time for Horn Schunck method = " << timeHS / count << "second(s)" << endl;
	cout << "Average time for Lucas Kanade method = " << timeLK / count << "second(s)" << endl;

	QString averageComputationTimeRecord = QString::fromUtf8(saveFolder.c_str()) + "\\file.csv";
	QFile file(averageComputationTimeRecord);
	if (file.open(QFile::WriteOnly | QFile::Truncate))
	{
		QTextStream stream(&file);
		stream << "Farneback,Liu,HS,LK\n"; // this writes first line with two columns
		stream << timeFarneback / count << "," << 
			timeLiu / count << "," <<
			timeHS / count << "," << 
			timeLK / count <<"\n"; // this writes first line with two columns
		file.close();
	}

	//system("pause");
	return EXIT_SUCCESS;
}